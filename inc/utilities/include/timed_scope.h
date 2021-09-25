#ifndef TIMED_SCOPE_H
#define TIMED_SCOPE_H

#include<iostream>

#ifndef WIN32
#include<sys/time.h>
#endif

namespace utilities {
  
    double get_time(){
#ifdef WIN32
        
#else
        struct timeval tv;
        gettimeofday( &tv, NULL );
        return double(tv.tv_sec) + 1e-6*double(tv.tv_usec);
#endif
    }
    
    class timed_scope {
    private:
        double  start_time;
    public:
        timed_scope(){
            start_time = get_time();
        }
        
        inline double timer() const {
            return get_time()-start_time;
        }
    };
    
};

std::ostream &operator<<( std::ostream &in, const utilities::timed_scope &ts ){
    in << ts.timer();
    return in;
}

#endif