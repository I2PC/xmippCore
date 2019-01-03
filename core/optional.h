/*
 * optional.h
 *
 *  Created on: Dec 12, 2018
 *      Author: david
 */

#ifndef CORE_OPTIONAL_H_
#define CORE_OPTIONAL_H_

#if __cplusplus >= 201703L
#include <optional>
namespace core::optional=std::optional
#else

#include <utility>
#include <assert.h>

namespace core {

template<typename T>
class optional {
public:
    optional() :
            p(new payload()) {
    }
    ;
    optional(const T &&t) :
            p(new payload_full(t)) {
    }
    ;

    optional(const T &t) :
            p(new payload_full(t)) {
    }
    ;

    ~optional() {
        delete p;
    }

    constexpr explicit operator bool() const {
        return p->has_value;
    }

    constexpr bool has_value() const {
        return p->has_value;
    }

    constexpr T value() const {
        assert(p->has_value);
        return static_cast<const payload_full*>(p)->t;
    }

private:

    struct payload {
        explicit payload(bool has_value = false) :
                has_value(has_value) {
        }
        ;
        const bool has_value;
    };

    struct payload_full: public payload {
        explicit payload_full(const T &t) :
                payload(true), t(std::move(t)) {
        }
        ;

        const T t;
    };

    const payload *p;

};
// optional

}// namespace

#endif 
#endif /* CORE_OPTIONAL_H_ */
