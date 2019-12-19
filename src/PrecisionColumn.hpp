//
// Created by Marc Suchard on 2019-12-10.
//

#ifndef ZIG_ZAG_PRECISION_COLUMN_HPP
#define ZIG_ZAG_PRECISION_COLUMN_HPP

namespace zz {

    class PrecisionColumnCallback {
    public:
        virtual const double *getColumn(int index) = 0;

        virtual void releaseColumn() = 0;
    };
}

#endif //ZIG_ZAG_PRECISION_COLUMN_HPP
