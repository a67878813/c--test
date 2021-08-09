#include <iostream>
#include <cmath>


int test1(void)
{
    int i = 0;
    int j = 0;
    std::cout << "begin" << std::endl;
    for (i = 1; i < 9; i++)
    {
        std::cout << i << std::endl;
    }
    return 0;
};

#define MACCCCCCCC 144

//123123s
// a
//appendix: Vector Operations
//==================
//------------------------------------------------------------------------
// Vector Class and vector functions
//------------------------------------------------------------------------
class Vector {
public:
    float x;
    float y;
    float z;

    Vector(void);
    Vector(float xi, float yi, float zi);

    float Magnitude(void);
    void  Normalize(void);
    void  Reverse(void);

    Vector& operator+=(Vector u);
    Vector& operator-=(Vector u);
    Vector& operator*=(float s);
    Vector& operator/=(float s);

    Vector operator-(void);

};

// Constructor
inline Vector::Vector(void)
{
    x = 0;
    y = 0;
    z = 0;
}

// Constructor
inline Vector::Vector(float xi, float yi, float zi)
{
    x = xi;
    y = yi;
    z = zi;
}
    
    
//====================================
inline    float Vector::Magnitude(void)
{
    return (float) sqrt(x*x + y*y + z*z);
}
    
float    const    tol = 0.0001f;
//====================================
inline    void  Vector::Normalize(void)
{
    float m = (float) sqrt(x*x + y*y + z*z);
    if(m <= tol) m = 1;
    x /= m;
    y /= m;
    z /= m;

    if (fabs(x) < tol) x = 0.0f;
    if (fabs(y) < tol) y = 0.0f;
    if (fabs(z) < tol) z = 0.0f;
}
    
    
//====================================

    
    
//====================================
inline    void  Vector::Reverse(void)
{
    x = -x;
    y = -y;
    z = -z;
}
    
    
//====================================
inline Vector& Vector::operator+=(Vector u)
{
    x += u.x;
    y += u.y;
    z += u.z;
    return *this;
}
    
    
//====================================
inline    Vector& Vector::operator-=(Vector u)
{
    x -= u.x;
    y -= u.y;
    z -= u.z;
    return *this;
}
    
    
//====================================
inline    Vector& Vector::operator*=(float s)
{
    x *= s;
    y *= s;
    z *= s;
    return *this;
}
    
    
//====================================
inline    Vector& Vector::operator/=(float s)
{
    x /= s;
    y /= s;
    z /= s;
    return *this;
}
    
    
//====================================
inline    Vector Vector::operator-(void)
{
    return Vector(-x, -y, -z);
}
    
    
//====================================
inline    Vector operator+(Vector u, Vector v)
{
    return Vector(u.x + v.x, u.y + v.y, u.z + v.z);
}
    
    
//====================================
inline    Vector operator-(Vector u, Vector v)
{
    return Vector(u.x - v.x, u.y - v.y, u.z - v.z);
}
    
    
//====================================
inline    Vector operator^(Vector u, Vector v)
{
    return Vector(    u.y*v.z - u.z*v.y,
                    -u.x*v.z + u.z*v.x,
                    u.x*v.y - u.y*v.x );
}
    
    
//====================================
// Vector dot product
inline    float operator*(Vector u, Vector v)
{
    return (u.x*v.x + u.y*v.y + u.z*v.z);
}
    
    
//====================================
inline    Vector operator*(float s, Vector u)
{
    return Vector(u.x*s, u.y*s, u.z*s);
}


inline    Vector operator*(Vector u, float s)
{
    return Vector(u.x*s, u.y*s, u.z*s);
}
    
    
//====================================
inline    Vector operator/(Vector u, float s)
{
    return Vector(u.x/s, u.y/s, u.z/s);
}
    
    
//====================================
inline    float TripleScalarProduct(Vector u, Vector v, Vector w)
{
    return float(    (u.x * (v.y*w.z - v.z*w.y)) +
                    (u.y * (-v.x*w.z + v.z*w.x)) +
                    (u.z * (v.x*w.y - v.y*w.x)) );
}
    
    





//123123123123123













// a end






typedef struct _PointMass{
    float mass;
    Vector designPosition;
    Vector correctedPosition;

} PointMass;

#define _NUMELEMENTS 10








int main()
{

    //test1();

    PointMass Elements[_NUMELEMENTS];



    int i ;
    float  TotalMass;
    Vector CombinedCG;
    Vector FirstMoment;




    for (i =0;i<_NUMELEMENTS;i++){
        Elements[i].mass = 1;
        Elements[i].designPosition.x = i ;
        std::cout <<" x= " << i << "\n";
    }











    TotalMass = 0;
    for(i=0 ;i<_NUMELEMENTS; i++){
        TotalMass += Elements[i].mass;
    };

    FirstMoment = Vector(0,0,0);
    for(i=0; i<_NUMELEMENTS; i++){
        FirstMoment += Elements[i].mass * Elements[i].designPosition;
    };
    // combined center of gravity
    CombinedCG = FirstMoment / TotalMass;
    std::cout <<" x= " << CombinedCG.x <<" y= " <<CombinedCG.y <<" z= " << CombinedCG.z << std::endl;


    
    std::cout <<"  "<< std::endl;

}
