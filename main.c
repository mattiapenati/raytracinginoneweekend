/*
 *            DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
 *                    Version 2, December 2004
 *
 * Copyright (C) 2004 Sam Hocevar <sam@hocevar.net>
 *
 * Everyone is permitted to copy and distribute verbatim or modified
 * copies of this license document, and changing it is allowed as long
 * as the name is changed.
 *
 *            DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
 *   TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION
 *
 *  0. You just DO WHAT THE FUCK YOU WANT TO.
 */

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#if 1
typedef float hFloat;
#define hSqrt sqrtf
#define hPow powf
#define hTan tanf
#else
typedef double hFloat;
#define hSqrt sqrt
#define hPow pow
#define hTan tan
#endif

typedef union hVec3f {
        struct { hFloat x, y, z; };
        struct { hFloat r, g ,b; };
} hVec3f;

void hVec3fZero(hVec3f *self);
void hVec3fInit(hVec3f *self, hFloat x, hFloat y, hFloat z);
void hVec3fAdd(hVec3f *self, const hVec3f *u, const hVec3f *v);
void hVec3fSub(hVec3f *self, const hVec3f *u, const hVec3f *v);
void hVec3fMul(hVec3f *self, const hVec3f *u, const hVec3f *v);
void hVec3fFma(hVec3f *self, hFloat a, const hVec3f *u, const hVec3f *v);
void hVec3fAdd1(hVec3f *self, hFloat a, const hVec3f *u);
void hVec3fScale(hVec3f *self, hFloat a, const hVec3f *u);
void hVec3fNorm(hVec3f *self, const hVec3f *other);
hFloat hVec3fDot(const hVec3f *u, const hVec3f *v);
hFloat hVec3fLength2(const hVec3f *u);
hFloat hVec3fLength(const hVec3f *u);
void hVec3fCross(hVec3f *self, const hVec3f *u, const hVec3f *v);
void hVec3fRandomSphere(hVec3f *self);


typedef struct hRay {
        hVec3f origin, direction;
} hRay;

void hRayEval(const hRay *self, hFloat t, hVec3f *ans);


typedef struct hCamera {
        hVec3f u, v, w;
        hVec3f origin, llcorner, hsize, vsize;
        hFloat radius;
} hCamera;

void hCameraInit(hCamera *self, const hVec3f *origin, const hVec3f *lookat,
                 const hVec3f *up, hFloat vfov, hFloat aspect, hFloat aperture,
                 hFloat distance);
void hCameraGenerateRay(const hCamera *self, hFloat u, hFloat v, hRay *ray);


typedef struct hIntersection {
        hFloat t;
        hVec3f point, normal;
        const void *material;
} hIntersection;


/* Shape interface */
typedef bool (*hShapeHitFn)(const void *, const hRay *, hFloat, hFloat,
                            hIntersection *);
typedef void (*hShapeClearFn)(void *);

typedef struct hShapeType {
        size_t size;

        hShapeHitFn hit;
        hShapeClearFn clear;
} hShapeType;

void * hShapeAlloc(const hShapeType *type);
void hShapeDestroy(void *self);
bool hShapeHit(const void *self, const hRay *ray, hFloat tmin, hFloat tmax,
               hIntersection *intersection);

void hShapeSetMaterial(void *self, const void *material);
const void * hShapeGetMaterial(const void *self);


/* Material interface */
typedef bool (*hMaterialScatterFn)(const void *, const hRay *,
                                   const hIntersection *, hVec3f *, hRay *);

typedef struct hMaterialType {
        size_t size;

        hMaterialScatterFn scatter;
} hMaterialType;

void * hMaterialAlloc(const hMaterialType *type);
void hMaterialDestroy(void *self);
bool hMaterialScatter(const void *self, const hRay *ray_in,
                      const hIntersection *intersection, hVec3f *attenuation,
                      hRay *ray_out);


/* Shapes */
typedef struct hSphere {
        hVec3f center;
        hFloat radius;
} hSphere;

void * hSphereNew(hFloat cx, hFloat cy, hFloat cz, hFloat radius);

bool hSphereHit(const hSphere *self, const hRay *ray, hFloat tmin, hFloat tmax,
                hIntersection *intersection);

static hShapeType hSphereType = {
        .size = sizeof(hSphere),
        .hit = (hShapeHitFn)hSphereHit,
};


typedef struct hShapeList {
        size_t size;
        void **items;
} hShapeList;

void * hShapeListNew(size_t size, void **items);

bool hShapeListHit(const hShapeList *self, const hRay *ray, hFloat tmin,
                   hFloat tmax, hIntersection *intersection);
void hShapeListClear(hShapeList *self);

static hShapeType hShapeListType = {
        .size = sizeof(hShapeList),
        .hit = (hShapeHitFn)hShapeListHit,
        .clear = (hShapeClearFn)hShapeListClear,
};


/* Materials */
typedef struct hLambertian {
        hVec3f albedo;
} hLambertian;

void * hLambertianNew(hFloat r, hFloat g, hFloat b);

bool hLambertianScatter(const hLambertian *self, const hRay *ray_in,
                        const hIntersection *intersection, hVec3f *attenuation,
                        hRay *ray_out);

static hMaterialType hLambertianType = {
        .size = sizeof(hLambertian),
        .scatter = (hMaterialScatterFn)hLambertianScatter,
};


typedef struct hMetal {
        hVec3f albedo;
        hFloat fuzz;
} hMetal;

void * hMetalNew(hFloat r, hFloat g, hFloat b, hFloat fuzz);

bool hMetalScatter(const hMetal *self, const hRay *ray_in,
                   const hIntersection *intersection, hVec3f *attenuation,
                   hRay *ray_out);

static hMaterialType hMetalType = {
        .size = sizeof(hMetal),
        .scatter = (hMaterialScatterFn)hMetalScatter,
};


typedef struct hDieletric {
        hFloat index;
} hDieletric;

void * hDieletricNew(hFloat index);

bool hDieletricScatter(const hDieletric *self, const hRay *ray_in,
                       const hIntersection *intersection, hVec3f *attenuation,
                       hRay *ray_out);

static hMaterialType hDieletricType = {
        .size = sizeof(hDieletric),
        .scatter = (hMaterialScatterFn)hDieletricScatter,
};


void
hPixelColor(hVec3f *color, void *world, hRay *ray_in, size_t depth)
{
        static hVec3f color0 = {.r = 1.0, .g = 1.0, .b = 1.0},
                      color1 = {.r = 0.5, .g = 0.7, .b = 1.0};

        hFloat t;
        hVec3f tmp, attenuation;
        hRay ray_out;
        hIntersection intersection;

        if (hShapeHit(world, ray_in, 1e-5, INFINITY, &intersection)) {
                if (depth < 50 && hMaterialScatter(intersection.material,
                                                   ray_in, &intersection,
                                                   &attenuation, &ray_out)) {
                        hPixelColor(&tmp, world, &ray_out, depth + 1);
                        hVec3fMul(color, &attenuation, &tmp);
                } else {
                        hVec3fZero(color);
                }
        } else {
                hVec3fNorm(&tmp, &ray_in->direction);
                t = (tmp.y + 1) / 2;
                hVec3fScale(color, 1 - t, &color0);
                hVec3fFma(color, t, &color1, color);
        }
}

void
hGammaCorrection(hVec3f *color, hFloat gamma)
{
        color->r = hPow(color->r, 1 / gamma);
        color->g = hPow(color->g, 1 / gamma);
        color->b = hPow(color->b, 1 / gamma);
}


#define NX 600
#define NY 300
unsigned char image[NY][NX][3];

int
main(int argc, char *argv[])
{
        size_t col, row, sample;
        float u, v;
        hVec3f color, tmp;
        hCamera camera;
        hRay ray;
        void *materials[4],
             *spheres[5],
             *world;

        /* camera */
        static const hVec3f lookfrom = {.x =3.0, .y = 3.0, .z = 2.0},
                            lookat = {.x = 0.0, .y = 0.0, .z = -1.0},
                            up = {.x = 0.0, .y = 1.0, .z = 0.0};
        static const hFloat aperture = 2.0,
                            distance = 5.1961524227;
        hCameraInit(&camera, &lookfrom, &lookat, &up, 20,
                    (hFloat)NX / (hFloat)NY, aperture, distance);
        static const size_t nsamples = 500;

        /* world */
        materials[0] = hLambertianNew(0.1, 0.2, 0.5);
        spheres[0] = hSphereNew(0, 0, -1, 0.5);
        hShapeSetMaterial(spheres[0], materials[0]);

        materials[1] = hLambertianNew(0.8, 0.8, 0.0);
        spheres[1] = hSphereNew(0, -100.5, -1, 100);
        hShapeSetMaterial(spheres[1], materials[1]);

        materials[2] = hMetalNew(0.8, 0.6, 0.2, 0.0);
        spheres[2] = hSphereNew(1, 0, -1, 0.5);
        hShapeSetMaterial(spheres[2], materials[2]);

        materials[3] = hDieletricNew(1.5);
        spheres[3] = hSphereNew(-1, 0, -1, 0.5);
        hShapeSetMaterial(spheres[3], materials[3]);

        spheres[4] = hSphereNew(-1, 0, -1, -0.45);
        hShapeSetMaterial(spheres[4], materials[3]);

        world = hShapeListNew(5, spheres);

        for (row = 0; row < NY; ++row)
                for (col = 0; col < NX; ++col) {
                        hVec3fZero(&color);
                        for (sample = 0; sample < nsamples; ++sample) {
                                u = (col + drand48()) / (NX - 1);
                                v = (row + drand48()) / (NY - 1);
                                hCameraGenerateRay(&camera, u, v, &ray);
                                hPixelColor(&tmp, world, &ray, 0);

                                hVec3fAdd(&color, &color, &tmp);
                        }
                        hVec3fScale(&color, 1. / nsamples, &color);
                        hGammaCorrection(&color, 2);

                        image[row][col][0] = 255 * color.r;
                        image[row][col][1] = 255 * color.g;
                        image[row][col][2] = 255 * color.b;
                }
        hShapeDestroy(world);
        hMaterialDestroy(materials[2]);
        hMaterialDestroy(materials[1]);
        hMaterialDestroy(materials[0]);


        printf("P3\n" "%d %d\n" "255\n", NX, NY);
        for (row = 0; row < NY; ++row)
                for (col = 0; col < NX; ++col)
                        printf("%d %d %d\n", image[NY - row - 1][col][0],
                                             image[NY - row - 1][col][1],
                                             image[NY - row - 1][col][2]);

        return EXIT_SUCCESS;
}


void
hVec3fZero(hVec3f *self)
{
        self->x = self->y = self->z = 0.0;
}

void
hVec3fInit(hVec3f *self, hFloat x, hFloat y, hFloat z)
{
        self->x = x;
        self->y = y;
        self->z = z;
}

void
hVec3fAdd(hVec3f *self, const hVec3f *u, const hVec3f *v)
{
        self->x = u->x + v->x;
        self->y = u->y + v->y;
        self->z = u->z + v->z;
}

void
hVec3fSub(hVec3f *self, const hVec3f *u, const hVec3f *v)
{
        self->x = u->x - v->x;
        self->y = u->y - v->y;
        self->z = u->z - v->z;
}

void
hVec3fMul(hVec3f *self, const hVec3f *u, const hVec3f *v)
{
        self->x = u->x * v->x;
        self->y = u->y * v->y;
        self->z = u->z * v->z;
}

void
hVec3fFma(hVec3f *self, hFloat a, const hVec3f *u, const hVec3f *v)
{
        self->x = a * u->x + v->x;
        self->y = a * u->y + v->y;
        self->z = a * u->z + v->z;
}

void
hVec3fAdd1(hVec3f *self, hFloat a, const hVec3f *u)
{
        self->x = a + u->x;
        self->y = a + u->y;
        self->z = a + u->z;
}

void
hVec3fScale(hVec3f *self, hFloat a, const hVec3f *u)
{
        self->x = u->x * a;
        self->y = u->y * a;
        self->z = u->z * a;
}

void
hVec3fNorm(hVec3f *self, const hVec3f *other)
{
        hVec3fScale(self, 1. / hVec3fLength(other), other);
}

hFloat
hVec3fDot(const hVec3f *u, const hVec3f *v)
{
        return (u->x * v->x) + (u->y * v->y) + (u->z * v->z);
}

hFloat
hVec3fLength2(const hVec3f *u)
{
        return hVec3fDot(u, u);
}

hFloat
hVec3fLength(const hVec3f *u)
{
        return hSqrt(hVec3fLength2(u));
}

void
hVec3fCross(hVec3f *self, const hVec3f *u, const hVec3f *v)
{
        self->x = u->y * v->z - u->z * v->y;
        self->y = u->z * v->x - u->x * v->z;
        self->z = u->x * v->y - u->y * v->x;
}

void
hVec3fRandomSphere(hVec3f *self)
{
        do {
                self->x = 2 * drand48() - 1;
                self->y = 2 * drand48() - 1;
                self->z = 2 * drand48() - 1;
        } while(hVec3fLength2(self) >= 1.0);
}


void
hRayEval(const hRay *self, hFloat t, hVec3f *ans)
{
        hVec3fFma(ans, t, &self->direction, &self->origin);
}


static void
hVec3fRandomUnitDisk(hVec3f *self)
{
        self->z = 0;
        do {
                self->x = 2 * drand48() - 1;
                self->y = 2 * drand48() - 1;
        } while(hVec3fLength2(self) >= 1.0);
}

void
hCameraInit(hCamera *self, const hVec3f *origin, const hVec3f *lookat,
            const hVec3f *up, hFloat vfov, hFloat aspect, hFloat aperture,
            hFloat distance)
{
        hFloat theta, hhsize, hvsize;

        hVec3fSub(&self->w, origin, lookat);

        hVec3fCross(&self->u, up, &self->w);
        hVec3fNorm(&self->u, &self->u);

        hVec3fCross(&self->v, &self->w, &self->u);
        hVec3fNorm(&self->v, &self->v);

        theta = vfov * M_PI / 180;
        hvsize = hTan(theta / 2);
        hhsize = aspect * hvsize;

        self->origin = *origin;

        hVec3fSub(&self->llcorner, origin, &self->w);
        hVec3fFma(&self->llcorner, -hhsize * distance, &self->u, &self->llcorner);
        hVec3fFma(&self->llcorner, -hvsize * distance, &self->v, &self->llcorner);

        hVec3fScale(&self->hsize, 2 * hhsize * distance, &self->u);
        hVec3fScale(&self->vsize, 2 * hvsize * distance, &self->v);

        self->radius = aperture / 2;
}

void
hCameraGenerateRay(const hCamera *self, hFloat u, hFloat v, hRay *ray)
{
        hVec3f rand;

        hVec3fRandomUnitDisk(&rand);
        rand.x *= self->radius * u;
        rand.y *= self->radius * v;

        hVec3fAdd(&ray->origin, &self->origin, &rand);

        hVec3fSub(&ray->direction, &self->llcorner, &ray->origin);
        hVec3fFma(&ray->direction, u, &self->hsize, &ray->direction);
        hVec3fFma(&ray->direction, v, &self->vsize, &ray->direction);
}


/* Shape interface */
typedef struct hShape {
        const hShapeType *type;
        const void *material;
        void *data;
} hShape;

static const hShape *
hShapeGetImpl(const void *self)
{
        return (const hShape *)((char *)self - offsetof(hShape, data));
}

void *
hShapeAlloc(const hShapeType *type)
{
        hShape *self = malloc(offsetof(hShape, data) + type->size);
        if (self) {
                self->type = type;
                self->material = NULL;
                return &(self->data);
        }
        return NULL;
}

void
hShapeDestroy(void *self)
{
        const hShape *shape = hShapeGetImpl(self);
        hShapeClearFn clear = shape->type->clear;
        if (clear)
                clear(self);

        free((void *)shape);
}

bool
hShapeHit(const void *self, const hRay *ray, hFloat tmin, hFloat tmax,
          hIntersection *intersection)
{
        const hShape *shape = hShapeGetImpl(self);
        hShapeHitFn hit = shape->type->hit;
        if (hit)
                return hit(self, ray, tmin, tmax, intersection);
        else
                return false;
}

void
hShapeSetMaterial(void *self, const void *material)
{
        hShape *shape = (hShape *)hShapeGetImpl(self);
        shape->material = material;
}

const void *
hShapeGetMaterial(const void *self)
{
        const hShape *shape = hShapeGetImpl(self);
        return shape->material;
}


/* Material interface */
typedef struct hMaterial {
        const hMaterialType *type;
        void *data;
} hMaterial;

static const hMaterial *
hMaterialGetImpl(const void *self)
{
        return (const hMaterial *)((char *)self - offsetof(hMaterial, data));
}

void *
hMaterialAlloc(const hMaterialType *type)
{
        hMaterial *self = malloc(offsetof(hMaterial, data) + type->size);
        if (self) {
                self->type = type;
                return &(self->data);
        }
        return NULL;
}

void
hMaterialDestroy(void *self)
{
        const hMaterial *material = hMaterialGetImpl(self);
        free((void *)material);
}

bool
hMaterialScatter(const void *self, const hRay *ray_in,
                 const hIntersection *intersection, hVec3f *attenuation,
                 hRay *ray_out)

{
        if (self) {
                const hMaterial *material = hMaterialGetImpl(self);
                hMaterialScatterFn scatter = material->type->scatter;
                return scatter(self, ray_in, intersection, attenuation, ray_out);
        }
        return false;
}

void * hMaterialAlloc(const hMaterialType *type);
void hMaterialDestroy(void *self);

/* Shapes */
void *
hSphereNew(hFloat cx, hFloat cy, hFloat cz, hFloat radius)
{
        hSphere *self = hShapeAlloc(&hSphereType);
        if (self) {
                self->center.x = cx;
                self->center.y = cy;
                self->center.z = cz;

                self->radius = radius;
        }
        return self;
}

bool
hSphereHit(const hSphere *self, const hRay *ray, hFloat tmin, hFloat tmax,
           hIntersection *intersection)
{
        hVec3f dir;
        hFloat a, b, c, disc, t;

        hVec3fSub(&dir, &ray->origin, &self->center);
        a = hVec3fDot(&ray->direction, &ray->direction);
        b = hVec3fDot(&dir, &ray->direction);
        c = hVec3fDot(&dir, &dir) - self->radius * self->radius;
        disc = b * b - a * c;

        if (disc > 0) {
                t = (-b - hSqrt(disc)) / a;
                if (tmin < t && t < tmax) {
                        intersection->t = t;
                        hRayEval(ray, t, &intersection->point);
                        hVec3fSub(&intersection->normal, &intersection->point,
                                  &self->center);
                        hVec3fScale(&intersection->normal, 1 / self->radius,
                                    &intersection->normal);
                        intersection->material = hShapeGetMaterial(self);
                        return true;
                }
                t = (-b + hSqrt(disc)) / a;
                if (tmin < t && t < tmax) {
                        intersection->t = t;
                        hRayEval(ray, t, &intersection->point);
                        hVec3fSub(&intersection->normal, &intersection->point,
                                  &self->center);
                        hVec3fScale(&intersection->normal, 1 / self->radius,
                                    &intersection->normal);
                        intersection->material = hShapeGetMaterial(self);
                        return true;
                }
        }
        return false;

}


void *
hShapeListNew(size_t size, void **items)
{
        hShapeList *self = hShapeAlloc(&hShapeListType);
        if (self) {
                self->size = size;
                self->items = malloc(size * sizeof(void *));
                memcpy(self->items, items, size * sizeof(void *));
        }
        return self;
}

bool
hShapeListHit(const hShapeList *self, const hRay *ray, hFloat tmin,
              hFloat tmax, hIntersection *intersection)
{
        bool ans;
        size_t idx;
        hFloat tclosest;

        ans = false;
        tclosest = tmax;
        for (idx = 0; idx < self->size; ++idx)
                if (hShapeHit(self->items[idx], ray, tmin, tclosest, intersection)) {
                        ans = true;
                        tclosest = intersection->t;
                }
        return ans;
}

void
hShapeListClear(hShapeList *self)
{
        size_t idx;

        for (idx = 0; idx < self->size; ++idx)
                hShapeDestroy(self->items[idx]);
        free(self->items);
}


/* Materials */
void *
hLambertianNew(hFloat r, hFloat g, hFloat b)
{
        hLambertian *self = hMaterialAlloc(&hLambertianType);
        if (self) {
                self->albedo.r = r;
                self->albedo.g = g;
                self->albedo.b = b;
        }
        return self;
}

bool
hLambertianScatter(const hLambertian *self, const hRay *ray_in,
                   const hIntersection *intersection, hVec3f *attenuation,
                   hRay *ray_out)
{
        hVec3f tmp;

        *attenuation = self->albedo;

        ray_out->origin = intersection->point;

        hVec3fRandomSphere(&tmp);
        hVec3fAdd(&ray_out->direction, &intersection->normal, &tmp);

        return true;
}


void *
hMetalNew(hFloat r, hFloat g, hFloat b, hFloat fuzz)
{
        hMetal *self = hMaterialAlloc(&hMetalType);
        if (self) {
                self->albedo.r = r;
                self->albedo.g = g;
                self->albedo.b = b;
                self->fuzz = fuzz;
        }
        return self;
}

static void
hVec3fReflect(hVec3f *self, const hVec3f *v, const hVec3f *n)
{
        hVec3fNorm(self, v);
        hVec3fFma(self, -2 * hVec3fDot(self, n), n, self);
}

bool
hMetalScatter(const hMetal *self, const hRay *ray_in,
              const hIntersection *intersection, hVec3f *attenuation,
              hRay *ray_out)
{
        hVec3f rand, reflected;

        *attenuation = self->albedo;

        hVec3fRandomSphere(&rand);
        hVec3fReflect(&reflected, &ray_in->direction, &intersection->normal);
        hVec3fFma(&reflected, self->fuzz, &rand, &reflected);

        ray_out->origin = intersection->point;
        ray_out->direction = reflected;

        return hVec3fDot(&ray_out->direction, &intersection->normal) > 0;
}


void *
hDieletricNew(hFloat index)
{
        hDieletric *self = hMaterialAlloc(&hDieletricType);
        if (self) {
                self->index = index;
        }
        return self;
}

static bool
hVec3fRefract(hVec3f *self, const hVec3f *v, const hVec3f *n, hFloat nu)
{
        hFloat dt, disc;

        hVec3fNorm(self, v);
        dt = hVec3fDot(self, n);
        disc = 1 - nu * nu * (1 - dt * dt);
        if (disc > 0) {
                hVec3fFma(self, -dt, n, self);
                hVec3fScale(self, nu, self);
                hVec3fFma(self, -hSqrt(disc), n, self);

                return true;
        }
        return false;
}

static hFloat
hSchlick(hFloat cosine, hFloat index)
{
        hFloat r0 = (1 - index) / (1 + index);
        r0 = r0 * r0;

        return r0 + (1 - r0) * hPow(1 - cosine, 5);
}

bool
hDieletricScatter(const hDieletric *self, const hRay *ray_in,
                  const hIntersection *intersection, hVec3f *attenuation,
                  hRay *ray_out)
{
        hVec3f reflected, refracted, normal;
        hFloat nu, cosine, prob;

        hVec3fInit(attenuation, 1.0, 1.0, 1.0);

        hVec3fReflect(&reflected, &ray_in->direction, &intersection->normal);

        normal = intersection->normal;
        nu = self->index;
        cosine = hVec3fDot(&ray_in->direction, &intersection->normal) /
                        hVec3fLength(&ray_in->direction);
        if (hVec3fDot(&ray_in->direction, &intersection->normal) > 0) {
                hVec3fScale(&normal, -1.0, &normal);
                cosine *= nu;
        } else {
                nu = 1.0 / nu;
                cosine *= -1;
        }

        if (hVec3fRefract(&refracted, &ray_in->direction, &normal, nu))
                prob = hSchlick(cosine, self->index);
        else
                prob = 1;

        ray_out->origin = intersection->point;
        if (drand48() <= prob)
                ray_out->direction = reflected;
        else
                ray_out->direction = refracted;

        return true;
}
