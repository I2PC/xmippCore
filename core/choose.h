#ifndef CORE_CHOOSE_H
#define CORE_CHOOSE_H

// TODO: write documentation & usage

template <bool flag, class typeTrue, class typeFalse>
struct choose;

template <class typeTrue, class typeFalse>
struct choose<true, typeTrue, typeFalse> {
   typedef typeTrue type;
};

template <class typeTrue, class typeFalse>
struct choose<false, typeTrue, typeFalse> {
   typedef typeFalse type;
};

#endif
