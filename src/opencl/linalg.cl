
/*!
 * This file provides structure and function definitions for the Vector2 and
 * Matrix2x2 types, which are vector and matrix types with fixed dimensions.
 * The operations defined for these types compute outputs directly without the
 * use of loops. These types are useful for any algorithm that operates on
 * pairwise data.
 *
 * Since OpenCL provides built-in vector types, Vector2 and Matrix2x2 are
 * defined in terms of these types. The following mapping is used to map
 * indices to xyzw:
 *
 *   ELEM(M, 0, 0) = M->x
 *   ELEM(M, 0, 1) = M->y
 *   ELEM(M, 1, 0) = M->z
 *   ELEM(M, 1, 1) = M->w
 */
typedef float2 Vector2;
typedef float4 Matrix2x2;






#define vectorInitZero(a) \
   (a)->x = 0; \
   (a)->y = 0;






#define vectorAdd(a, b) \
   (a)->x += (b)->x; \
   (a)->y += (b)->y;






#define vectorAddScaled(a, c, b) \
   (a)->x += (c) * (b)->x; \
   (a)->y += (c) * (b)->y;






#define vectorSubtract(a, b) \
   (a)->x -= (b)->x; \
   (a)->y -= (b)->y;






#define vectorScale(a, c) \
   (a)->x *= (c); \
   (a)->y *= (c);






#define vectorDot(a, b) \
   ((a)->x * (b)->x + (a)->y * (b)->y)






#define SQR(x) ((x)*(x))
#define vectorDiffNorm(a, b) \
   sqrt(SQR((a)->x - (b)->x) + SQR((a)->y - (b)->y))






#define matrixInitIdentity(M) \
   (M)->x = 1; \
   (M)->y = 0; \
   (M)->z = 0; \
   (M)->w = 1;






#define matrixInitZero(M) \
   (M)->x = 0; \
   (M)->y = 0; \
   (M)->z = 0; \
   (M)->w = 0;






#define matrixScale(A, c) \
   (A)->x *= (c); \
   (A)->y *= (c); \
   (A)->z *= (c); \
   (A)->w *= (c);






#define matrixInverse(A, B, det) \
   *det = (A)->x * (A)->w - (A)->y * (A)->z; \
   (B)->x = +(A)->w / (*det); \
   (B)->y = -(A)->y / (*det); \
   (B)->z = -(A)->z / (*det); \
   (B)->w = +(A)->x / (*det);






#define matrixProduct(A, x_, b) \
   (b)->x = (A)->x * (x_)->x + (A)->y * (x_)->y; \
   (b)->y = (A)->z * (x_)->x + (A)->w * (x_)->y;






#define matrixAddOuterProduct(A, c, x_) \
   (A)->x += (c) * (x_)->x * (x_)->x; \
   (A)->y += (c) * (x_)->x * (x_)->y; \
   (A)->z += (c) * (x_)->y * (x_)->x; \
   (A)->w += (c) * (x_)->y * (x_)->y;
