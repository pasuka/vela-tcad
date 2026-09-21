"""Generate branch-preserving IALMob scalar/explicit-partial kernels (SymPy 1.14).

Differentiate a named scalar expression graph offline. CSE is local to each
branch: never hoist singular expressions out of zero/absent-scattering guards.
The seven SI directions and the stationary-screening envelope convention match
IalMobilityEvaluation.h. Runtime code has no overloaded AD arithmetic.
"""
import argparse
import re
from pathlib import Path
import sympy as s

FIELDS='mu3 mu2 transitionScale coulombDamping damping phonon3d phononNumerator phononTemperaturePower roughnessExponent roughnessNumerator'.split()


class Graph:
    def __init__(self, differentiated):
        self.diff=differentiated;self.lines=[];self.der={};self.env={};self.serial=0
    def sym(self,name):
        if name not in self.env:self.env[name]=s.Symbol(name,real=True)
        return self.env[name]
    def parse(self,expr):
        for name in re.findall(r'\b[a-zA-Z_]\w*\b',expr):
            if name not in ('exp','log','sqrt'):self.sym(name)
        return s.sympify(expr,locals=dict(self.env,exp=s.exp,log=s.log,sqrt=s.sqrt))
    def code(self,expr):
        return re.sub(r'\bp_(\w+)\b',r'p.\1',s.ccode(expr))
    def load(self,name,value,derivatives):
        self.sym(name);self.lines.append('    const Real %s=%s;'%(name,value))
        self.der[name]=[s.sympify(x) for x in derivatives] if self.diff else [s.S.Zero]*7
    def expression(self,name,expr,declare=True,guard=None):
        expr=self.parse(expr) if isinstance(expr,str) else expr
        deps=[x for x in sorted(expr.free_symbols,key=str) if str(x) in self.der]
        ds=[sum((s.diff(expr,x)*self.der[str(x)][k] for x in deps),s.S.Zero) for k in range(7)] if self.diff else []
        # Keep primal arithmetic identical for residual and Jacobian calls.
        # Derivative CSE must not rewrite the value or hoist cusp slopes.
        groups=([expr],ds)
        reduced=[]
        for i,group in enumerate(groups):
            if not group:continue
            pairs,values=s.cse(group,s.numbered_symbols('c%d_'%self.serial),order='canonical');self.serial+=1
            if guard and i==1:
                # Partial common expressions live only inside the selected branch.
                for k in range(7):self.lines.append('    Real %s_d%d=0.;'%(name,k))
                self.lines.append('    if('+guard+') {')
            for x,v in pairs:self.lines.append('    const Real %s=%s;'%(x,self.code(v)))
            if guard and i==1:
                for k,v in enumerate(values):self.lines.append('    %s_d%d=%s;'%(name,k,self.code(v)))
                self.lines.append('    }')
            else:reduced.extend(values)
        self.lines.append('    '+('Real ' if declare else '')+name+'='+self.code(reduced[0])+';')
        self.sym(name)
        if self.diff:
            if not guard:
                for k,v in enumerate(reduced[1:]):self.lines.append('    '+('Real ' if declare else '')+'%s_d%d=%s;'%(name,k,self.code(v)))
            self.der[name]=[self.sym(name+'_d%d'%k) for k in range(7)]
        else:self.der[name]=[s.S.Zero]*7
        return name
    def raw(self,line):self.lines.append('    '+line)
    def result(self,name):
        return 'generatedDual('+name+', {'+', '.join(map(self.code,self.der[name]))+'})' if self.diff else name
    def power(self,name,base,exponent):
        return self.expression(name,base+'**('+exponent+')',guard=base+'>0.')


def preparation(d):
    g=Graph(d)
    for i,name in enumerate(('ndSI','naSI','nSI','hSI','fieldSI','distanceSI','temperature')):
        g.load(name,'state[%d]'%i,[int(i==k) for k in range(7)])
    for name,expr in [('nd','ndSI/1e6'),('na','naSI/1e6'),('n','nSI/1e6'),('h','hSI/1e6'),('distance','distanceSI*100'),('t','temperature/300')]:g.expression(name,expr)
    for name,a,b in [('carrier','n','h'),('other','h','n'),('inv','na','nd'),('acc','nd','na')]:
        g.load(name,'electron?%s:%s'%(a,b),[0]*7)
        if d:
            g.der[name]=[s.Symbol(name+'_d%d'%k) for k in range(7)]
            for k in range(7):g.raw('const Real %s_d%d=electron?%s:%s;'%(name,k,g.code(g.der[a][k]),g.code(g.der[b][k])))
    for name,doping,ref,coef in [('ndStar','nd','nRefD','cRefD'),('naStar','na','nRefA','cRefA')]:
        g.expression(name,doping+'*(1+('+doping+'/p_'+ref+')**2/(1+p_'+coef+'*('+doping+'/p_'+ref+')**2))')
    g.expression('nsc','ndStar+naStar+other');g.expression('mu3','INFINITY');g.expression('screenG','1')
    g.raw('if(nsc>0.) {')
    g.expression('screening','t*t/((2.459/3.97e13)*nsc**(2./3.)+(3.828/(p_mass*1.36e20))*(carrier+other))')
    g.load('clamped','screening>pMin?screening:pMin',[0]*7)
    if d:
        g.der['clamped']=[s.Symbol('clamped_d%d'%k) for k in range(7)]
        for k in range(7):g.raw('const Real clamped_d%d=screening>pMin?screening_d%d:0.;'%(k,k))
    g.expression('screenG','1-.89233/(.41372+clamped*(t/p_mass)**.28227)**.19778+.005978/(clamped*(p_mass/t)**.72169)**1.80618',False)
    g.expression('pPower','screening**.6478')
    g.expression('f','(.7643*pPower+2.2999+6.5502*p_mass/p_otherMass)/(pPower+2.3670-.8552*p_mass/p_otherMass)')
    g.load('effBase','electron?ndStar:naStar',[0]*7);g.load('effOther','electron?naStar:ndStar',[0]*7)
    if d:
        for name,a,b in [('effBase','ndStar','naStar'),('effOther','naStar','ndStar')]:
            g.der[name]=[s.Symbol(name+'_d%d'%k) for k in range(7)]
            for k in range(7):g.raw('const Real %s_d%d=electron?%s:%s;'%(name,k,g.code(g.der[a][k]),g.code(g.der[b][k])))
    g.expression('effective','effBase+screenG*effOther+other/f')
    g.expression('mu3','(p_muMax*p_muMax/(p_muMax-p_muMin))*t**(3*p_alpha-1.5)*nsc/effective*(p_nRef/nsc)**p_alpha+(p_muMax*p_muMin/(p_muMax-p_muMin))*t**(-.5)*(carrier+other)/effective',False)
    g.raw('}')
    for name,doping,suffix in [('muInv','inv','Inv'),('muAcc','acc','Acc')]:
        g.expression(name,'INFINITY');g.raw('if('+doping+'>0.) {')
        g.expression(name+'Norm','carrier/p_nScRef');g.power(name+'Carrier',name+'Norm','p_nu0'+suffix)
        g.expression(name+'A','p_d1'+suffix+'*t**p_alpha1'+suffix+'*'+name+'Carrier/('+doping+'/p_nDopRef)**p_nu1'+suffix)
        g.expression(name+'B','p_d2'+suffix+'*t**p_alpha2'+suffix+'/('+doping+'/p_nDopRef)**p_nu2'+suffix)
        g.load(name+'Scale','std::max(std::abs('+name+'A),std::abs('+name+'B))',[0]*7)
        g.raw('if('+name+'Scale==0.) {');g.expression(name,'0',False);g.raw('} else {')
        g.expression(name,name+'Scale*sqrt(('+name+'A/'+name+'Scale)**2+('+name+'B/'+name+'Scale)**2)',False)
        g.raw('}');g.raw('}')
    g.expression('accG','INFINITY');g.raw('if(std::isfinite(muAcc)) {');g.expression('accG','muAcc*screenG',False);g.raw('}')
    for name,x in [('inverseInv','muInv'),('inverseAcc','accG')]:
        g.expression(name,'0');g.raw('if(!std::isinf('+x+')) {');g.expression(name,'1/'+x,False);g.raw('}')
    g.expression('inverse2','inverseInv+inverseAcc');g.expression('mu2','INFINITY')
    g.raw('if(inverse2!=0.) {');g.expression('mu2','1/inverse2',False);g.raw('}')
    for name,expr in [('transitionScale','p_S/temperature'),('coulombDamping','exp(-distance/p_lCritC)'),('damping','exp(-distance/p_lCrit)'),('phonon3d','p_muMax*t**(-p_theta)'),('phononNumerator','p_C*(na+nd+p_N2)**p_lambda'),('phononTemperaturePower','t**p_k'),('roughnessExponent','p_A+p_alphaSr*(n+h)/(na+nd+p_N1)**p_nu'),('roughnessNumerator','(na+nd+p_N2)**p_lambdaSr')]:
        # N2=0 is allowed. Match the reference selected power slope at zero,
        # including exponent=0; do not evaluate 0*pow(0,-1).
        g.expression(name,expr,guard='na+nd+p.N2>0.' if name in ('phononNumerator','roughnessNumerator') else None)
    g.raw('return {'+', '.join(g.result(x) for x in FIELDS)+'};')
    return '\n'.join(g.lines)


def evaluation(d):
    g=Graph(d)
    for name in FIELDS:g.load(name,'q.'+name+('.value' if d else ''),[s.Symbol('q.'+name+'.derivative[%d]'%k) for k in range(7)])
    g.load('fieldSI','field_SI',[int(k==4) for k in range(7)])
    g.expression('field','fieldSI/100');g.power('fieldPower','field','2./3.')
    g.expression('transition','transitionScale*fieldPower-p_transitionP');g.expression('f','0')
    g.raw('if(transition<=700.) {');g.expression('f','1/(1+exp(transition))',False);g.raw('}')
    g.expression('weight','coulombDamping*(1-f)');g.expression('muC','INFINITY')
    g.raw('if(weight==0.) {');g.expression('muC','mu3',False);g.raw('} else if(weight==1.) {');g.expression('muC','mu2',False)
    g.raw('} else if(std::isfinite(mu3)&&std::isfinite(mu2)) {');g.expression('muC','(1-weight)*mu3+weight*mu2',False);g.raw('}')
    g.expression('muPh','phonon3d');g.expression('muSr','INFINITY');g.raw('if(field>0. && damping>0.) {')
    g.expression('muPh2','p_B/field+phononNumerator/(field**(1./3.)*phononTemperaturePower)')
    g.expression('muPh','1/(damping/muPh2+1/phonon3d)',False)
    g.expression('inverseSr','field**roughnessExponent/p_delta+field**3/p_eta')
    g.expression('muSr','roughnessNumerator/inverseSr',False);g.raw('}')
    for name,x in [('invC','muC'),('invPh','muPh'),('invSr','muSr')]:
        g.expression(name,'0');g.raw('if(!std::isinf('+x+')) {');g.expression(name,'1/'+x,False);g.raw('}')
    g.expression('mobility','1/(invC+invPh+damping*invSr)')
    outputs=[]
    for name in ('mobility','mu3','mu2','muC','muPh','muSr'):
        g.expression('out_'+name,name+'*1e-4');outputs.append('out_'+name)
    g.raw('return {'+', '.join(g.result(x) for x in outputs)+'};')
    return '\n'.join(g.lines)


def generate():
    parts=['// Generated by scripts/generate_ialmob_symbolic.py; SymPy '+s.__version__+'. DO NOT EDIT.',
        '#pragma once','#include "vela/physics/detail/IalMobilityEvaluation.h"',
        'namespace vela::ial_detail {','using std::pow; using std::exp; using std::log; using std::sqrt;',
        'inline Dual generatedDual(Real value,std::array<Real,7> d) { Dual out(value);out.derivative=d;return out; }']
    for d in (False,True):
        typ='Dual' if d else 'Real';suffix='Partials' if d else 'Value'
        parts += ['inline Preparation<'+typ+'> generatedPrepare'+suffix+'(const std::array<Real,7>& state,const IalMobilityParameters& p,bool electron,Real pMin) {',preparation(d),'}',
            'inline std::array<'+typ+',6> generatedEvaluate'+suffix+'(const Preparation<'+typ+'>& q,Real field_SI,const IalMobilityParameters& p) {',evaluation(d),'}']
    parts += ['} // namespace vela::ial_detail','']
    return '\n'.join(parts)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('include/vela/physics/detail/IalMobilityGenerated.h'))
    parser.add_argument('--check',action='store_true');args=parser.parse_args();text=generate()
    if args.check:
        if args.output.read_text(encoding='utf-8')!=text:raise SystemExit('Generated kernel differs')
    else:args.output.write_text(text,encoding='utf-8')
