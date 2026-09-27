// vela: native launcher for the sdevice-compatible front end (contract C9).
//
//   vela device <args...>  ->  <python> <bin>/vela_sdevice.pyz device <args...>
//
// The child runs in the same working directory with the same standard
// handles; every argument is forwarded unchanged and the child's exit code is
// returned as is.  <bin> is <exe dir>/bin when that holds vela_sdevice.pyz,
// otherwise the launcher's own directory.  Environment overrides:
//   VELA_SDEVICE_PYZ  path of the orchestrator zipapp
//   VELA_PYTHON       interpreter (path or bare command name)
// Interpreter order: VELA_PYTHON -> <bin>/python/ (bundled) -> python3/python
// on PATH (warning).  Standard library only; no third-party dependencies.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cwctype>
#include <filesystem>
#include <initializer_list>
#include <optional>
#include <string>
#include <system_error>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
#include <cerrno>
#include <unistd.h>
#endif

#ifndef VELA_VERSION
#define VELA_VERSION "0.0.0"
#endif
#ifndef VELA_GIT_DESCRIBE
#define VELA_GIT_DESCRIBE ""
#endif

namespace {

namespace fs = std::filesystem;

#ifdef _WIN32
using NativeChar = wchar_t;
using NativeString = std::wstring;
#define VELA_T(literal) L##literal
constexpr NativeChar kPathListSeparator = L';';
#else
using NativeChar = char;
using NativeString = std::string;
#define VELA_T(literal) literal
constexpr NativeChar kPathListSeparator = ':';
#endif

constexpr const char* kSdeviceSyntax = "T-2022.03";

constexpr const char* kUsage =
    "Usage: vela <command> [arguments...]\n"
    "       vela --version | -v\n"
    "       vela --help | -h\n"
    "\n"
    "Commands:\n"
    "  device    Run a Sentaurus Device (sdevice) command file with Vela.\n"
    "            Forwards to: <python> <bin>/vela_sdevice.pyz device [arguments...]\n"
    "\n"
    "Environment:\n"
    "  VELA_PYTHON        Python interpreter to use (default: <bin>/python/, then PATH)\n"
    "  VELA_SDEVICE_PYZ   Orchestrator archive (default: <bin>/vela_sdevice.pyz)\n";

std::string toUtf8(const NativeString& text)
{
#ifdef _WIN32
    if (text.empty()) {
        return {};
    }
    const int size = WideCharToMultiByte(CP_UTF8, 0, text.data(), static_cast<int>(text.size()),
                                         nullptr, 0, nullptr, nullptr);
    if (size <= 0) {
        return "?";
    }
    std::string out(static_cast<std::size_t>(size), '\0');
    WideCharToMultiByte(CP_UTF8, 0, text.data(), static_cast<int>(text.size()), out.data(), size,
                        nullptr, nullptr);
    return out;
#else
    return text;
#endif
}

std::string toUtf8(const fs::path& path)
{
    return toUtf8(path.native());
}

void printErr(const std::string& message)
{
    std::fputs(message.c_str(), stderr);
    std::fflush(stderr);
}

std::optional<NativeString> getEnv(const NativeChar* name)
{
#ifdef _WIN32
    const DWORD size = GetEnvironmentVariableW(name, nullptr, 0);
    if (size == 0) {
        return std::nullopt;
    }
    std::wstring value(size, L'\0');
    const DWORD written = GetEnvironmentVariableW(name, value.data(), size);
    if (written == 0 || written >= size) {
        return std::nullopt;
    }
    value.resize(written);
    return value;
#else
    const char* value = std::getenv(name);
    if (value == nullptr || *value == '\0') {
        return std::nullopt;
    }
    return NativeString(value);
#endif
}

void setEnv(const NativeChar* name, const NativeString& value)
{
#ifdef _WIN32
    SetEnvironmentVariableW(name, value.c_str());
#else
    setenv(name, value.c_str(), 1);
#endif
}

void unsetEnv(const NativeChar* name)
{
#ifdef _WIN32
    SetEnvironmentVariableW(name, nullptr);
#else
    unsetenv(name);
#endif
}

bool isRegularFile(const fs::path& path)
{
    std::error_code ec;
    return fs::is_regular_file(path, ec);
}

bool isExecutableFile(const fs::path& path)
{
#ifdef _WIN32
    return isRegularFile(path);
#else
    return isRegularFile(path) && access(path.c_str(), X_OK) == 0;
#endif
}

#ifdef _WIN32
bool containsIgnoreCase(const std::wstring& haystack, const std::wstring& needle)
{
    std::wstring lowered(haystack);
    for (wchar_t& ch : lowered) {
        ch = static_cast<wchar_t>(std::towlower(static_cast<wint_t>(ch)));
    }
    return lowered.find(needle) != std::wstring::npos;
}
#endif

fs::path executablePath(const NativeChar* argv0)
{
    fs::path found;
#ifdef _WIN32
    std::wstring buffer(MAX_PATH, L'\0');
    for (;;) {
        const DWORD length =
            GetModuleFileNameW(nullptr, buffer.data(), static_cast<DWORD>(buffer.size()));
        if (length == 0) {
            break;
        }
        if (length < buffer.size()) {
            buffer.resize(length);
            found = fs::path(buffer);
            break;
        }
        buffer.resize(buffer.size() * 2);
    }
#else
    std::error_code linkError;
    const fs::path link = fs::read_symlink("/proc/self/exe", linkError);
    if (!linkError && !link.empty()) {
        found = link;
    }
#endif
    if (found.empty()) {
        std::error_code absError;
        found = fs::absolute(fs::path(argv0), absError);
        if (absError) {
            found = fs::path(argv0);
        }
    }
    std::error_code canonicalError;
    const fs::path canonical = fs::weakly_canonical(found, canonicalError);
    return canonicalError ? found : canonical;
}

struct Layout {
    fs::path exeDir;
    fs::path binDir;
    fs::path pyz;  // empty when not found
    std::vector<fs::path> searched;
};

Layout discoverLayout(const fs::path& exePath)
{
    Layout layout;
    layout.exeDir = exePath.parent_path();
    const fs::path nestedBin = layout.exeDir / VELA_T("bin");
    std::error_code ec;
    layout.binDir = fs::is_directory(nestedBin, ec) ? nestedBin : layout.exeDir;

    if (const auto overridePath = getEnv(VELA_T("VELA_SDEVICE_PYZ"))) {
        const fs::path requested(*overridePath);
        layout.searched.push_back(requested);
        if (isRegularFile(requested)) {
            layout.pyz = requested;
        }
        return layout;
    }
    for (const fs::path& candidate : {nestedBin / VELA_T("vela_sdevice.pyz"),
                                      layout.exeDir / VELA_T("vela_sdevice.pyz")}) {
        layout.searched.push_back(candidate);
        if (isRegularFile(candidate)) {
            layout.pyz = candidate;
            layout.binDir = candidate.parent_path();
            return layout;
        }
    }
    return layout;
}

std::vector<fs::path> pathEntries()
{
    std::vector<fs::path> entries;
    const auto value = getEnv(VELA_T("PATH"));
    if (!value) {
        return entries;
    }
    NativeString::size_type start = 0;
    while (start <= value->size()) {
        NativeString::size_type end = value->find(kPathListSeparator, start);
        if (end == NativeString::npos) {
            end = value->size();
        }
        NativeString entry = value->substr(start, end - start);
        if (entry.size() >= 2 && entry.front() == VELA_T('"') && entry.back() == VELA_T('"')) {
            entry = entry.substr(1, entry.size() - 2);
        }
        if (!entry.empty()) {
            entries.emplace_back(entry);
        }
        start = end + 1;
    }
    return entries;
}

std::optional<fs::path> findOnPath(const std::vector<NativeString>& names)
{
    const std::vector<fs::path> directories = pathEntries();
    for (const NativeString& name : names) {
        for (const fs::path& directory : directories) {
#ifdef _WIN32
            // Store app-execution aliases are not runnable interpreters.
            if (containsIgnoreCase(directory.native(), L"\\microsoft\\windowsapps")) {
                continue;
            }
#endif
            const fs::path candidate = directory / name;
            if (isExecutableFile(candidate)) {
                return candidate;
            }
        }
    }
    return std::nullopt;
}

std::vector<fs::path> bundledInterpreterCandidates(const fs::path& binDir)
{
    const fs::path root = binDir / VELA_T("python");
#ifdef _WIN32
    return {root / L"python.exe"};
#else
    return {root / "bin" / "python3", root / "bin" / "python", root / "python3", root / "python"};
#endif
}

struct Interpreter {
    fs::path path;
    bool bundled = false;
    bool fromPath = false;
};

std::optional<Interpreter> resolveInterpreter(const Layout& layout)
{
    if (const auto value = getEnv(VELA_T("VELA_PYTHON"))) {
        const fs::path requested(*value);
        if (isExecutableFile(requested)) {
            return Interpreter{requested, false, false};
        }
        if (!requested.has_parent_path()) {
            std::vector<NativeString> names{requested.native()};
#ifdef _WIN32
            if (requested.extension().empty()) {
                names.push_back(requested.native() + L".exe");
            }
#endif
            if (const auto found = findOnPath(names)) {
                return Interpreter{*found, false, false};
            }
        }
        printErr("vela: VELA_PYTHON='" + toUtf8(requested) +
                 "' is not an executable file; ignoring it\n");
    }
    for (const fs::path& candidate : bundledInterpreterCandidates(layout.binDir)) {
        if (isExecutableFile(candidate)) {
            return Interpreter{candidate, true, false};
        }
    }
#ifdef _WIN32
    const std::vector<NativeString> names{L"python3.exe", L"python.exe"};
#else
    const std::vector<NativeString> names{"python3", "python"};
#endif
    if (const auto found = findOnPath(names)) {
        return Interpreter{*found, false, true};
    }
    return std::nullopt;
}

void prepareChildEnvironment(const Layout& layout, const Interpreter& interpreter)
{
    setEnv(VELA_T("PYTHONNOUSERSITE"), VELA_T("1"));
    if (interpreter.bundled) {
        unsetEnv(VELA_T("PYTHONPATH"));
        unsetEnv(VELA_T("PYTHONHOME"));
    }
    NativeString path = layout.binDir.native();
    if (const auto existing = getEnv(VELA_T("PATH"))) {
        path += kPathListSeparator;
        path += *existing;
    }
    setEnv(VELA_T("PATH"), path);
}

#ifdef _WIN32
// Standard CommandLineToArgvW / CRT quoting: quote when needed, double the
// backslashes that precede a quote or the closing quote, escape quotes.
void appendQuoted(std::wstring& commandLine, const std::wstring& argument)
{
    if (!argument.empty() && argument.find_first_of(L" \t\n\v\"") == std::wstring::npos) {
        commandLine += argument;
        return;
    }
    commandLine += L'"';
    std::size_t backslashes = 0;
    for (const wchar_t ch : argument) {
        if (ch == L'\\') {
            ++backslashes;
            continue;
        }
        if (ch == L'"') {
            commandLine.append(backslashes * 2 + 1, L'\\');
        } else {
            commandLine.append(backslashes, L'\\');
        }
        backslashes = 0;
        commandLine += ch;
    }
    commandLine.append(backslashes * 2, L'\\');
    commandLine += L'"';
}

int runChild(const fs::path& interpreter, const std::vector<NativeString>& arguments)
{
    std::wstring commandLine;
    appendQuoted(commandLine, interpreter.native());
    for (const std::wstring& argument : arguments) {
        commandLine += L' ';
        appendQuoted(commandLine, argument);
    }

    STARTUPINFOW startup{};
    startup.cb = sizeof(startup);
    startup.dwFlags = STARTF_USESTDHANDLES;
    startup.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
    startup.hStdOutput = GetStdHandle(STD_OUTPUT_HANDLE);
    startup.hStdError = GetStdHandle(STD_ERROR_HANDLE);
    PROCESS_INFORMATION process{};
    if (!CreateProcessW(interpreter.c_str(), commandLine.data(), nullptr, nullptr, TRUE, 0,
                        nullptr, nullptr, &startup, &process)) {
        printErr("vela: E202 cannot start '" + toUtf8(interpreter) + "' (Windows error " +
                 std::to_string(GetLastError()) + ")\n");
        return 2;
    }
    CloseHandle(process.hThread);
    // Ctrl+C is handled by the child; the launcher only reports its exit code.
    // Must follow CreateProcessW: the ignore flag would otherwise be inherited.
    SetConsoleCtrlHandler(nullptr, TRUE);
    WaitForSingleObject(process.hProcess, INFINITE);
    DWORD exitCode = 2;
    if (!GetExitCodeProcess(process.hProcess, &exitCode)) {
        exitCode = 2;
    }
    CloseHandle(process.hProcess);
    return static_cast<int>(exitCode);
}
#else
int runChild(const fs::path& interpreter, const std::vector<NativeString>& arguments)
{
    std::string program = interpreter.native();
    std::vector<std::string> storage(arguments);
    std::vector<char*> argv;
    argv.reserve(storage.size() + 2);
    argv.push_back(program.data());
    for (std::string& argument : storage) {
        argv.push_back(argument.data());
    }
    argv.push_back(nullptr);
    std::fflush(nullptr);
    execv(program.c_str(), argv.data());
    printErr("vela: E202 cannot execute '" + program + "': " + std::strerror(errno) + "\n");
    return 2;
}
#endif

int printVersion()
{
    const char* describe = VELA_GIT_DESCRIBE;
    std::printf("vela %s (%s) sdevice-syntax %s\n", VELA_VERSION,
                (*describe != '\0') ? describe : "unknown", kSdeviceSyntax);
    std::fflush(stdout);
    return 0;
}

int launcherMain(int argc, NativeChar** argv)
{
    if (argc < 2) {
        printErr("vela: E201 no command given\n");
        printErr(kUsage);
        return 2;
    }
    const NativeString command = argv[1];
    if (command == VELA_T("--version") || command == VELA_T("-v")) {
        return printVersion();
    }
    if (command == VELA_T("--help") || command == VELA_T("-h")) {
        std::fputs(kUsage, stdout);
        std::fflush(stdout);
        return 0;
    }
    if (command != VELA_T("device")) {
        printErr("vela: E202 unknown command '" + toUtf8(command) + "'\n");
        printErr(kUsage);
        return 2;
    }

    const Layout layout = discoverLayout(executablePath(argv[0]));
    if (layout.pyz.empty()) {
        std::string message = "vela: E202 vela_sdevice.pyz not found; searched:";
        for (const fs::path& candidate : layout.searched) {
            message += "\n  " + toUtf8(candidate);
        }
        message += "\nSet VELA_SDEVICE_PYZ to the orchestrator archive.\n";
        printErr(message);
        return 2;
    }

    const std::optional<Interpreter> interpreter = resolveInterpreter(layout);
    if (!interpreter) {
        printErr("vela: E202 no Python interpreter available "
                 "(set VELA_PYTHON, install <bin>/python/, or add python3 to PATH)\n");
        return 2;
    }
    if (interpreter->fromPath) {
        printErr("vela: using Python from PATH (" + toUtf8(interpreter->path) + ")\n");
    }

    prepareChildEnvironment(layout, *interpreter);

    std::vector<NativeString> arguments;
    arguments.reserve(static_cast<std::size_t>(argc));
    arguments.push_back(layout.pyz.native());
    for (int i = 1; i < argc; ++i) {
        arguments.emplace_back(argv[i]);
    }
    return runChild(interpreter->path, arguments);
}

}  // namespace

#ifdef _WIN32
int wmain(int argc, wchar_t** argv)
{
    return launcherMain(argc, argv);
}
#else
int main(int argc, char** argv)
{
    return launcherMain(argc, argv);
}
#endif
