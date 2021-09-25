#ifndef FLIP_VEC2_H
#define FLIP_VEC2_H

#include<cmath>
#include<iostream>

namespace geom {
        
    template< typename real >
    class vec2 {
    private:
        union {
            struct {
                real m_x;
                real m_y;
            };
            real m_v[2];
        };
    public:
        vec2( ){
            set( 0.0, 0.0 );
        }
        
        vec2( const real x, const real y ){
            set( x, y );
        }
        
        inline const real *ptr() const {
            return m_v;
        }
        
        inline void set( real x, real y ){
            m_x = x; m_y = y;
        }
        
        inline real dot( const vec2<real> &in ) const {
            return m_x*in.m_x + m_y*in.m_y;
        }
        
        inline real cross( const vec2<real> &in ) const {
            return m_x*in.m_y - m_y*in.m_x;
        }
        
        inline real length_squared() const {
            return m_x*m_x + m_y*m_y;
        }
        
        inline real length() const {
            return sqrt( m_x*m_x + m_y*m_y );
        }
        
        inline real normalize(){
            real len = length();
            m_x /= len; m_y /= len;
            return len;
        }
        
        inline vec2<real> component_parallel_to( const vec2<real> &b ) const {
            real blen2 = b.length_squared();
            return b*dot(b)/blen2;
        }
        
        inline vec2<real> component_perpendicular_to( const vec2<real> &b ) const {
            return operator-( component_parallel_to( b ) );
        }
        
        inline real &x(){
            return m_x;
        }
        
        inline const real &x() const {
            return m_x;
        }
        
        inline real &y(){
            return m_y;
        }
        
        inline const real &y() const {
            return m_y;
        }
        
        inline real &operator[]( const int id ){
            return m_v[id];
        }
        
        inline const real &operator[]( const int id ) const {
            return m_v[id];
        }
        
        inline vec2<real> operator-() const {
            return vec2<real>( -m_x, -m_y );
        }
        
        inline vec2<real> operator+( const vec2<real> b ) const {
            return vec2<real>( m_x+b.m_x, m_y+b.m_y );
        }
        
        inline vec2<real> operator-( const vec2<real> b ) const {
            return vec2<real>( m_x-b.m_x, m_y-b.m_y );
        }
        
        inline vec2<real> operator*( const real b ) const {
            return vec2<real>( m_x*b, m_y*b );
        }
        
        inline vec2<real> operator/( const real b ) const {
            return vec2<real>( m_x/b, m_y/b );
        }
        
        inline vec2<real> operator+=( const vec2<real> b ){
            vec2<real> tmp( m_x+b.m_x, m_y+b.m_y );
            *this=tmp;
            return tmp;
        }
        
        inline vec2<real> operator-=( const vec2<real> b ){
            vec2<real> tmp( m_x-b.m_x, m_y-b.m_y );
            *this=tmp;
            return tmp;
        }
        
        inline vec2<real> operator*=( const real b ){
            vec2<real> tmp( m_x*b, m_y*b );
            *this=tmp;
            return tmp;
        }
        
        inline vec2<real> operator/=( const real b ){
            vec2<real> tmp( m_x/b, m_y/b );
            *this=tmp;
            return tmp;
        }
    };
    
    template<typename real>
    vec2<real> operator*( const float &b, const vec2<real> &a ){
        return vec2<real>( b*a.x(), b*a.y() );
    }
    
    template<typename real>
    vec2<real> operator*( const double &b, const vec2<real> &a ){
        return vec2<real>( b*a.x(), b*a.y() );
    }

};

template<typename real>
std::ostream& operator<<( std::ostream &in, const geom::vec2<real> &val ){
    in << "[" << val.x() << ", " << val.y() << "]";
    return in;
}

#endif