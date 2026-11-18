#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include <filesystem>
#include <string>

namespace {
std::wstring quote(const std::filesystem::path& path)
{
    return L"\"" + path.wstring() + L"\"";
}

std::wstring environment(const wchar_t* name)
{
    const DWORD length = GetEnvironmentVariableW(name, nullptr, 0);
    std::wstring value(length, L'\0');
    if (length) {
        value.resize(GetEnvironmentVariableW(name, value.data(), length));
    }
    return value;
}
}

int WINAPI wWinMain(HINSTANCE, HINSTANCE, PWSTR, int)
{
    std::wstring program(32768, L'\0');
    const DWORD length = GetModuleFileNameW(nullptr, program.data(), static_cast<DWORD>(program.size()));
    if (!length || length >= program.size()) {
        return 1;
    }
    program.resize(length);
    const auto root = std::filesystem::path(program).parent_path();
    const auto pythonDir = root / L"python" / L"bin";
    const auto python = pythonDir / L"pythonw.exe";
    const auto runtime = root / L"vela";
    const auto mkl = runtime / L"bin";

    SetEnvironmentVariableW(L"PATH", (pythonDir.wstring() + L";" + mkl.wstring()
                                    + L";" + environment(L"PATH")).c_str());
    SetEnvironmentVariableW(L"PYTHONHOME", (root / L"python").c_str());
    SetEnvironmentVariableW(L"PYTHONPATH", nullptr);
    SetEnvironmentVariableW(L"PYTHONNOUSERSITE", L"1");
    SetEnvironmentVariableW(L"PYTHONDONTWRITEBYTECODE", L"1");
    SetEnvironmentVariableW(L"VELA_PYTHON", (pythonDir / L"python.exe").c_str());
    SetEnvironmentVariableW(L"VELA_BIN", mkl.c_str());
    SetEnvironmentVariableW(L"VELA_MKL_RT_LIBRARY", (mkl / L"mkl_rt.3.dll").c_str());

    auto command = quote(python) + L" " + quote(root / L"app" / L"vela_workbench.py");
    STARTUPINFOW startup{};
    startup.cb = sizeof(startup);
    PROCESS_INFORMATION process{};
    if (!CreateProcessW(python.c_str(), command.data(), nullptr, nullptr, FALSE, 0,
                        nullptr, root.c_str(), &startup, &process)) {
        MessageBoxW(nullptr, L"Could not start the bundled Python runtime. Extract the entire ZIP first.",
                    L"Semiconductor Device Simulator", MB_OK | MB_ICONERROR);
        return 1;
    }
    CloseHandle(process.hThread);
    CloseHandle(process.hProcess);
    return 0;
}
