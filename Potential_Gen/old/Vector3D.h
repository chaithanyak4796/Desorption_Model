/***************************************************************************************************
****************************************************************************************************

    Vector3D.h - defines a class for vectors with x, y, and z components

****************************************************************************************************
***************************************************************************************************/
class Vector3D
{
    public:
        Vector3D(); // constructor
        double x; // x component
        double y; // y component
        double z; // z component
};

Vector3D::Vector3D()
{
    x = 0;
    y = 0;
    z = 0;
}
