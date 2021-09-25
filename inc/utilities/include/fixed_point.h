#ifndef FIXED_POINT_H
#define FIXED_POINT_H
 
/**
 @file fixed_point.h
 @author James Gregson (james.gregson@gmail.com)
 @brief A simple templated fixed-point arithmetic class
*/

#include<iostream>

namespace utilities {
		
	/**
	 @brief fixed point number class, takes three template parameters. These paramters are:
		- DECBITS     the number of fractional bits to use
		- STOREDTYPE  the underlying integer type of the number
		- COMPTYPE    the integer type to be used for computations
	*/
	template< uint16_t DECBITS=10, typename STOREDTYPE=int64_t, typename COMPTYPE=int64_t >
	class fixed {
	public:
		typedef fixed<DECBITS,STOREDTYPE,COMPTYPE> fixed_type;
	private:
		STOREDTYPE	m_value;
	public:
		fixed( ) : m_value(0) { }
		fixed( const float a  ) : m_value( a*float(COMPTYPE(1)<<DECBITS) ){ }
		fixed( const double a ) : m_value( a*double(COMPTYPE(1)<<DECBITS) ) { }
		fixed( const fixed_type &in ) : m_value(in.m_value) { }
		
		inline operator int64_t() const { return m_value>>DECBITS; }
		inline operator double() const { return double(m_value)/double(COMPTYPE(1)<<DECBITS); }
		inline operator float() const { return double(m_value)/double(COMPTYPE(1)<<DECBITS); }
		
		inline double to_double() const { return double(m_value)/double(COMPTYPE(1)<<DECBITS); }
		
		inline std::ostream& print( std::ostream &os ) const {
			os << this->to_double();
			return os;
		}
				
		inline fixed_type operator+=( const fixed_type& in ){ m_value+=in.m_value; return *this; }
		inline fixed_type operator-=( const fixed_type& in ){ m_value-=in.m_value; return *this; }
		inline fixed_type operator*=( const fixed_type& in ){
			m_value = STOREDTYPE((COMPTYPE(m_value)*COMPTYPE(in.m_value))>>DECBITS);
			return *this;
		}
		inline fixed_type operator/=( const fixed_type& in ){
			m_value = STOREDTYPE((COMPTYPE(m_value)<<COMPTYPE(DECBITS))/COMPTYPE(in.m_value));
			return *this;
		}
		
		inline fixed_type operator-() const { fixed_type b; b.m_value=-m_value; return b; }
		
		inline friend fixed_type operator+( fixed_type a, const fixed_type& b ){ return a+=b; }
		inline friend fixed_type operator-( fixed_type a, const fixed_type& b ){ return a-=b; }
		inline friend fixed_type operator*( fixed_type a, const fixed_type& b ){ return a*=b; }
		inline friend fixed_type operator/( fixed_type a, const fixed_type& b ){ return a/=b; }
		
		inline bool operator==( const fixed_type& in ) const { return m_value==in.m_value; }
		inline bool operator!=( const fixed_type& in ) const { return m_value!=in.m_value; }
		inline bool operator<(  const fixed_type& in ) const { return m_value<in.m_value; }
		inline bool operator>(  const fixed_type& in ) const { return m_value>in.m_value; }
		inline bool operator<=( const fixed_type& in ) const { return m_value<=in.m_value; }
		inline bool operator>=( const fixed_type& in ) const { return m_value>=in.m_value; }
	};
	
	template< uint16_t DECBITS, typename STORED, typename COMPTYPE >
	std::istream& operator>>( std::istream &is, utilities::fixed<DECBITS,STORED,COMPTYPE> &in ){
		double v;
		is >> v;
		in = utilities::fixed<DECBITS,STORED,COMPTYPE>(v);
		return is;
	}
			
	template< uint16_t DECBITS, typename STORED, typename COMPTYPE >
	std::ostream& operator<<( std::ostream &os, const utilities::fixed<DECBITS,STORED,COMPTYPE> &in ){
		return in.print(os);
	}

		
};

		
#endif