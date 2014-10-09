#ifndef MATHTOOL_H
#define MATHTOOL_H

template<class E>
const E & min(const E & a, const E & b)
{
    return (a<b)?a:b;
}

template<class E>
const E & max(const E & a, const E & b)
{
    return (a>b)?a:b;
}

#endif
