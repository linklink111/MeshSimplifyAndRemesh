#ifndef FLIP_VEC2_H
#define FLIP_VEC2_H

#include<cmath>
#include<algorithm>
#include<iostream>
#include<drand48.h>

namespace geom {
    
    template< typename real >
    class vec3 {
    private:
		real m_v[3];
    public:
		vec3(){
			set(real(0.0),real(0.0),real(0.0));
		}
		
		vec3( const vec3<real> &in ){
			x() = in.x(); y()=in.y(); z()=in.z();
		}
        
        vec3( const real _x, const real _y=0.0, const real _z=0.0 ){
            set( _x, _y, _z );
        }
        
        inline const real *ptr() const {
            return m_v;
        }
        
        inline void set( const real _x, const real _y, const real _z ){
            x() = _x; y() = _y; z() = _z;
        }
        
        inline real dot( const vec3<real> &in ) const {
            return x()*in.x() + y()*in.y() + z()*in.z();
        }
        
        inline vec3<real> cross( const vec3<real> &in ) const {
            return vec3<real>( y()*in.z()-z()*in.y(), z()*in.x()-x()*in.z(), x()*in.y()-y()*in.x() );
        }
        
        inline real length_squared() const {
            return this->dot(*this);
        }
        
        inline real length() const {
            return sqrt( this->dot(*this) );
        }
        
        inline real normalize( const real safe_tol=real(0.0) ){
            real len = std::max( safe_tol, length() );
            x() /= len; y() /= len; z() /= len;
            return len;
        }
        
        inline vec3<real> component_parallel_to( const vec3<real> &b, const real safe_tol=real(0.0) ) const {
            real blen2 = std::max( safe_tol, b.length_squared() );
            return b*dot(b)/blen2;
        }
        
        inline vec3<real> component_perpendicular_to( const vec3<real> &b ) const {
            return operator-( component_parallel_to( b ) );
        }
        
        inline real &x(){
            return m_v[0];
        }
        
        inline const real &x() const {
            return m_v[0];
        }
        
        inline real &y(){
            return m_v[1];
        }
        
        inline const real &y() const {
            return m_v[1];
        }
        
        inline real &z(){
            return m_v[2];
        }
        
        inline const real &z() const {
            return m_v[2];
        }
        
        inline real &operator[]( const int id ){
            return m_v[id];
        }
        
        inline const real &operator[]( const int id ) const {
            return m_v[id];
        }
        
        inline vec3<real> operator-() const {
            return vec3<real>( -x(), -y(), -z() );
        }
        
        inline vec3<real> operator+( const vec3<real> b ) const {
            return vec3<real>( x()+b.x(), y()+b.y(), z()+b.z() );
        }
        
        inline vec3<real> operator-( const vec3<real> b ) const {
            return vec3<real>( x()-b.x(), y()-b.y(), z()-b.z() );
        }
        
        inline vec3<real> operator*( const real b ) const {
            return vec3<real>( x()*b, y()*b, z()*b );
        }
        
        inline vec3<real> operator/( const real b ) const {
            return vec3<real>( x()/b, y()/b, z()/b );
        }
        
        inline vec3<real> operator+=( const vec3<real> b ){
            vec3<real> tmp( x()+b.x(), y()+b.y(), z()+b.z() );
            *this=tmp;
            return tmp;
        }
        
        inline vec3<real> operator-=( const vec3<real> b ){
            vec3<real> tmp( x()-b.x(), y()-b.y(), z()-b.z() );
            *this=tmp;
            return tmp;
        }
        
        inline vec3<real> operator*=( const real b ){
            vec3<real> tmp( x()*b, y()*b, z()*b );
            *this=tmp;
            return tmp;
        }
        
        inline vec3<real> operator/=( const real b ){
            vec3<real> tmp( x()/b, y()/b, z()/b );
            *this=tmp;
            return tmp;
        }
        
        inline vec3<real> min( const vec3<real> &in ) const {
            return vec3<real>( std::min( x(), in[0] ), std::min( y(), in[1] ), std::min( z(), in[2] ) );
        }
        
        inline vec3<real> max( const vec3<real> &in ) const {
            return vec3<real>( std::max( x(), in[0] ), std::max( y(), in[1] ), std::max( z(), in[2] ) );
        }
        
        inline static vec3<real> random( const real size=1.0 ){
            return vec3<real>( size*(drand48()-0.5), size*(drand48()-0.5), size*(drand48()-0.5) );
        }

    };
    
    template<typename real>
    vec3<real> operator*( const float &b, const vec3<real> &a_ ){
        return vec3<real>( b*a_.x(), b*a_.y(), b*a_.z() );
    }
    
    template<typename real>
    vec3<real> operator*( const double &b, const vec3<real> &a_ ){
        return vec3<real>( b*a_.x(), b*a_.y(), b*a_.z() );
    }
	
	template<typename real>
	std::ostream& operator<<( std::ostream &in, const vec3<real> &val ){
		in << "[" << val.x() << "," << val.y() << "," << val.z() << "]";
		return in;
	}

    
};


#endif