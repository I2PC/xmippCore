#ifndef CORE_MEMORY_H
#define CORE_MEMORY_H

#include <memory>

namespace MemHelpers {

// This unit contains memory helpers, mainly rewrites of <memory> file to use
// newer constructions in older C++.

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

}//end namespace MemHelpers

#endif
