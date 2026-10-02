#include <catch2/catch_test_macros.hpp>
#include "simplemos_linear_refinement.hpp"
#include <limits>

using namespace simplemos_diagnostic;

TEST_CASE("High precision dot residual retains a weak term lost in double summation") {
    vela::SparseMatrixd J(1,2);J.insert(0,0)=1e20;J.insert(0,1)=-1e20;
    vela::VectorXd F(1),x(2);F<<1.;x<<1.,1.;
    auto d=highPrecisionDefect(J,F,x);
    REQUIRE(d.values(0)==1.);
    REQUIRE(d.backwardMax>0.);
    REQUIRE(d.backwardMax<1e-19);
}

TEST_CASE("Refinement recovers a missing weak equation update with original row weights") {
    vela::SparseMatrixd J(3,3);J.insert(0,0)=2;J.insert(1,1)=1e-30;J.insert(2,2)=3;
    vela::VectorXd F(3),step(3),w(3);F<<-2.,-2.5e-47,-6.;step<<1.,0.,2.;w<<1.,1e12,1.;
    vela::LinearSolver solver;int observations=0;
    auto result=refine(J,F,w,step,solver,4,[&](int k,const auto&,const Defect& d){
        ++observations;if(k==0)REQUIRE(d.values(1)==F(1));
    });
    REQUIRE(observations==5);
    REQUIRE(std::abs(result(1)/2.5e-17-1)<1e-14);
    REQUIRE(result(0)==1.);REQUIRE(result(2)==2.);
    REQUIRE(std::abs(highPrecisionDefect(J,F,result).values(1)/F(1))<1e-14);
}

TEST_CASE("An already exact solve is invariant under fixed refinement") {
    vela::SparseMatrixd J(2,2);J.insert(0,0)=4.;J.insert(0,1)=-2.;J.insert(1,0)=1.;J.insert(1,1)=2.;
    vela::VectorXd F(2),step(2),w(2);step<<.5,.25;F<<-1.5,-1.;w<<8.,.125;
    vela::LinearSolver solver;
    auto result=refine(J,F,w,step,solver,4,[](int,const auto&,const Defect& d){REQUIRE(d.values.norm()==0.);});
    REQUIRE((result-step).norm()==0.);
}

TEST_CASE("Zero corrections preserves the supplied raw direction") {
    vela::SparseMatrixd J(2,2);J.setIdentity();
    vela::VectorXd F=vela::VectorXd::Ones(2),step=vela::VectorXd::Zero(2),w=F;
    vela::LinearSolver solver;int count=0;
    auto result=refine(J,F,w,step,solver,0,[&](int,const auto&,const auto&){++count;});
    REQUIRE(count==1);REQUIRE(result.norm()==0.);
}
