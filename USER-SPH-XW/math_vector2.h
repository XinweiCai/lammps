#ifndef MATH_VECTOR_H_
#define MATH_VECTOR_H_

#include <array>
#include <cmath>
#include <cassert>
#include <type_traits>
#include <algorithm>

template<typename SCALAR, uint D>
struct vector {
/*---------------------------------------------------------------------------
								 Container
---------------------------------------------------------------------------*/
protected:
	SCALAR x[D];
public:
	using TYPE_ = SCALAR;
	static const uint D_ = D;

	constexpr inline uint d() const { return D_; }
	constexpr inline uint size() const { return D_; }

	// default constructor
	inline vector() {}
	// proxy for constructing from scalar OR pointer
	template<typename T> inline vector(T const initializer) {
		init< T, std::is_pointer<T>::value >( *this, initializer );
	}
protected:
	// construct from scalar constant
	template<typename T, bool is_ptr = false> struct init {
		inline init( vector &v, T const scalar ) {
			for(uint i = 0 ; i < D ; i++) v[i] = scalar;
		}
	};
	// construct from pointer to array
	template<typename P> struct init<P,true> {
		inline init( vector &v, P const ptr ) {
			for(uint i = 0 ; i < D ; i++) v[i] = ptr[i];
		}
	};
public:
	// assign-construct
	inline vector( vector const &other ) {
		for(uint i = 0; i < D; i++) x[i] = other.x[i];
	}
	// construct from parameter pack
	// 'head' differentiate it from constructing from vector expression
	// 'SCALAR const head' does NOT work because it confuses with init from single scalar (with tail... = void)
	template<typename H, typename ...T> inline vector(H const head, T const ... tail ) {
		std::array<TYPE_,D_> s( { static_cast<SCALAR>(head), static_cast<SCALAR>(tail)... } );
		for(uint i = 0 ; i < D ; i++) x[i] = s[i];
	}
	// copy-assign
	inline vector & operator =( vector const &u ) {
		for(uint i = 0 ; i < D ; i++) x[i] = u[i];
		return *this;
	}
	// TODO: move construct
	// TODO: move assign

	// TODO: is this really a good way?
	template<typename T> inline operator vector<T,D> () {
		vector<T,D> v;
		for(uint i=0;i<D;i++) v[i] = static_cast<T>( x[i] );
		return v;
	}
	template<typename T> inline vector( vector<T,D> const &other ) {
		for(uint i = 0; i < D; i++) x[i] = static_cast<T>( other[i] );
	}

	// point must be assignable, while other expressions may not
	inline SCALAR      & operator [] (uint i)       { return x[i]; }
	inline SCALAR const& operator [] (uint i) const { return x[i]; }

	// STL-style direct data accessor
	inline SCALAR      * data()       { return x; }
	inline SCALAR const* data() const { return x; }

/*---------------------------------------------------------------------------
						 Operator Overloads
---------------------------------------------------------------------------*/
	// OP-Assign operators
	inline vector & operator += ( vector const &u ) {
		for(uint i = 0 ; i < D ; i++) x[i] += u[i];
		return *this;
	}
	inline vector & operator -= ( vector const &u ) {
		for(uint i = 0 ; i < D ; i++) x[i] -= u[i];
		return *this;
	}
	inline vector & operator *= ( vector const &u ) {
		for(uint i = 0 ; i < D ; i++) x[i] *= u[i];
		return *this;
	}
	inline vector & operator /= ( vector const &u ) {
		for(uint i = 0 ; i < D ; i++) x[i] /= u[i];
		return *this;
	}
	// Vector-Scalar operators
	inline vector & operator += ( SCALAR const u ) {
		for(uint i = 0 ; i < D ; i++) x[i] += u;
		return *this;
	}
	inline vector & operator -= ( SCALAR const u ) {
		for(uint i = 0 ; i < D ; i++) x[i] -= u;
		return *this;
	}
	inline vector & operator *= ( SCALAR const u ) {
		for(uint i = 0 ; i < D ; i++) x[i] *= u;
		return *this;
	}
	inline vector & operator /= ( SCALAR const u ) {
		for(uint i = 0 ; i < D ; i++) x[i] /= u;
		return *this;
	}

	friend inline vector operator + ( vector const &u, vector const &v ) {
		vector w;
		for(uint i=0;i<D;i++) w[i] = u[i] + v[i];
		return w;
	}
	friend inline vector operator - ( vector const &u, vector const &v ) {
		vector w;
		for(uint i=0;i<D;i++) w[i] = u[i] - v[i];
		return w;
	}
	friend inline vector operator - ( vector const &u ) {
		vector w;
		for(uint i=0;i<D;i++) w[i] = -u[i];
		return w;
	}
	friend inline vector operator + ( vector const &u, SCALAR a ) {
		vector w;
		for(uint i=0;i<D;i++) w[i] = u[i] + a;
		return w;
	}
	friend inline vector operator - ( vector const &u, SCALAR a ) {
		vector w;
		for(uint i=0;i<D;i++) w[i] = u[i] - a;
		return w;
	}
	friend inline vector operator * ( vector const &u, vector const &v ) {
		vector w;
		for(uint i=0;i<D;i++) w[i] = u[i] * v[i];
		return w;
	}
	friend inline vector operator * ( vector const &u, SCALAR a ) {
		vector w;
		for(uint i=0;i<D;i++) w[i] = u[i] * a;
		return w;
	}
	friend inline vector operator * ( SCALAR a, vector const &u ) {
		return u * a;
	}
	friend inline vector operator / ( vector const &u, vector const &v ) {
		vector w;
		for(uint i=0;i<D;i++) w[i] = u[i] / v[i];
		return w;
	}
	friend inline vector operator / ( vector const &u, SCALAR a ) {
		vector w;
		for(uint i=0;i<D;i++) w[i] = u[i] / a;
		return w;
	}
	friend inline vector operator / ( SCALAR a, vector const &u ) {
		vector w;
		for(uint i=0;i<D;i++) w[i] = a / u[i];
		return w;
	}
	template<typename OP> friend inline vector apply( vector const &u, OP const &o ) {
		vector w;
		for(uint i=0;i<D;i++) w[i] = o( u[i] );
		return w;
	}
	template<typename OP> friend inline vector apply( vector const &u, vector const &v, OP const &o ) {
		vector w;
		for(uint i=0;i<D;i++) w[i] = o( u[i], v[i] );
		return w;
	}

/*---------------------------------------------------------------------------
						 Math functions
---------------------------------------------------------------------------*/
	// generic reduction template
	template<class OP> friend inline SCALAR reduce( vector const &u, OP const & op ) {
		SCALAR core( u[0] );
		for(uint i = 1 ; i < D ; i++) core = op( core, u[i] );
		return core;
	}
	// biggest element within a vector
	friend inline SCALAR max( vector const &u ) {
		return reduce( u, [](SCALAR a, SCALAR b){return a>b?a:b;} );
	}
	// smallest element within a vector
	friend inline SCALAR min( vector const &u ) {
		return reduce( u, [](SCALAR a, SCALAR b){return a<b?a:b;} );
	}
	// sum of elements
	friend inline SCALAR sum( vector const &u ) {
		return reduce( u, [](SCALAR a, SCALAR b){return a+b;} );
	}
	// mean of elements
	friend inline SCALAR mean( vector const &u ) {
		return sum(u) / double(D);
	}
	// inner product
	friend inline SCALAR dot( vector const &u, vector const &v ) {
		return sum( u * v );
	}
	// square of L2 norm
	friend inline SCALAR normsq( vector const &u ) {
		return sum( u * u );
	}
	// L2 norm
	friend inline SCALAR norm( vector const &u ) {
		return std::sqrt( normsq(u) );
	}
};

template<typename SCALAR> inline
vector<SCALAR,3U> cross( vector<SCALAR,3U> const &u, vector<SCALAR,3U> const &v ) {
	vector<SCALAR,3U> x;
	for(uint i=0;i<3;i++)  x[i] = u[(i+1U)%3U] * v[(i+2U)%3U] - u[(i+2U)%3U] * v[(i+1U)%3U];
	return x;
}

namespace functors {
template<typename T> struct min { inline T operator () ( T const &u, T const &v ) const { return std::min<T>( u, v ); } };
template<typename T> struct max { inline T operator () ( T const &u, T const &v ) const { return std::max<T>( u, v ); } };
}

#endif
