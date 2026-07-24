#include <filesystem>
#include <stdexcept>

static std::filesystem::path findProjectRootFromExe(const char* argv0)
{
    namespace fs = std::filesystem;
    fs::path exeDir = fs::absolute(argv0).parent_path();

    // 从 exeDir 往上找，直到找到含有 cases 文件夹的目录
    for (fs::path p = exeDir; !p.empty(); p = p.parent_path())
    {
        if (fs::is_directory(p / "cases"))
            return p;
        // 防止死循环
        if (p == p.root_path()) break;
    }
    throw std::runtime_error("Cannot locate project root (cases folder not found upward from exe path).");
}
