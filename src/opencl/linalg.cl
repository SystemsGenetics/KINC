
typedef union
{
   float s[2];
   float2 v2;
} Vector2;

typedef union
{
   float s[4];
   float4 v4;
} Matrix2x2;






#define ELEM(M, i, j) ((M)->s[(i) * 2 + (j)])






#define vectorInitZero(a) \
   (a)->s[0] = 0; \
   (a)->s[1] = 0;






#define vectorAdd(a, b) \
   (a)->s[0] += (b)->s[0]; \
   (a)->s[1] += (b)->s[1];






#define vectorAddScaled(a, c, b) \
   (a)->s[0] += (c) * (b)->s[0]; \
   (a)->s[1] += (c) * (b)->s[1];






#define vectorSubtract(a, b) \
   (a)->s[0] -= (b)->s[0]; \
   (a)->s[1] -= (b)->s[1];






#define vectorScale(a, c) \
   (a)->s[0] *= (c); \
   (a)->s[1] *= (c);






#define vectorDot(a, b) \
   ((a)->s[0] * (b)->s[0] + (a)->s[1] * (b)->s[1])






#define SQR(x) ((x)*(x))
#define vectorDiffNorm(a, b) \
   sqrt(SQR((a)->s[0] - (b)->s[0]) + SQR((a)->s[1] - (b)->s[1]))






#define matrixInitIdentity(M) \
   ELEM(M, 0, 0) = 1; \
   ELEM(M, 0, 1) = 0; \
   ELEM(M, 1, 0) = 0; \
   ELEM(M, 1, 1) = 1;






#define matrixInitZero(M) \
   ELEM(M, 0, 0) = 0; \
   ELEM(M, 0, 1) = 0; \
   ELEM(M, 1, 0) = 0; \
   ELEM(M, 1, 1) = 0;






#define matrixAddScaled(A, c, B) \
   ELEM(A, 0, 0) += (c) * ELEM(B, 0, 0); \
   ELEM(A, 0, 1) += (c) * ELEM(B, 0, 1); \
   ELEM(A, 1, 0) += (c) * ELEM(B, 1, 0); \
   ELEM(A, 1, 1) += (c) * ELEM(B, 1, 1);






#define matrixScale(A, c) \
   ELEM(A, 0, 0) *= (c); \
   ELEM(A, 0, 1) *= (c); \
   ELEM(A, 1, 0) *= (c); \
   ELEM(A, 1, 1) *= (c);






#define matrixInverse(A, B, det) \
   *det = ELEM(A, 0, 0) * ELEM(A, 1, 1) - ELEM(A, 0, 1) * ELEM(A, 1, 0); \
   ELEM(B, 0, 0) = +ELEM(A, 1, 1) / (*det); \
   ELEM(B, 0, 1) = -ELEM(A, 0, 1) / (*det); \
   ELEM(B, 1, 0) = -ELEM(A, 1, 0) / (*det); \
   ELEM(B, 1, 1) = +ELEM(A, 0, 0) / (*det);






#define matrixProduct(A, x, b) \
   (b)->s[0] = ELEM(A, 0, 0) * (x)->s[0] + ELEM(A, 0, 1) * (x)->s[1]; \
   (b)->s[1] = ELEM(A, 1, 0) * (x)->s[0] + ELEM(A, 1, 1) * (x)->s[1];






#define matrixOuterProduct(a, b, C) \
   ELEM(C, 0, 0) = (a)->s[0] * (b)->s[0]; \
   ELEM(C, 0, 1) = (a)->s[0] * (b)->s[1]; \
   ELEM(C, 1, 0) = (a)->s[1] * (b)->s[0]; \
   ELEM(C, 1, 1) = (a)->s[1] * (b)->s[1];
