#include <bits/stdc++.h>
#ifndef M_PI
#define M_PI acos(-1)
#endif
constexpr double eps = 1e-4; 
struct Vec
{
    Vec(double x_ = 0, double y_ = 0, double z_ = 0){
        x = x_;y = y_;z = z_;
    }
    double x,y,z; // position or color (r,g,b)
    Vec operator+(const Vec &b) const { return {x + b.x, y + b.y, z + b.z}; }
    Vec operator-(const Vec &b) const { return {x - b.x, y - b.y, z - b.z}; }
    Vec operator*(double b) const { return {x * b, y * b, z * b}; }
    Vec mult(const Vec &b) const { return {x * b.x, y * b.y, z * b.z}; }      // 计算颜色用的
    Vec &norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); } // 当前向量单位化
    double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }    // 内积
    Vec operator%(const Vec &b) const
    { // 求this和b的叉乘
        return {y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x};
    }
    double greatest() const // 获取RGB最大值
    {
        if ((x > y) && (x > z)) return x;
        if (y > z) return y;
        return z;
    }
};

struct Ray
{
    Vec o;
    Vec d; // 光线定义是原点o，以及单位方向向量d，射线就由o+td表达，t是一个参数，表达了到o点的距离
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};

enum class Refl
{
    DIFF,
    SPEC,
    REFR
}; // material types, used in radiance()，依次漫反射、镜面反射、折射

class Object
{
public:
    Object()
    {
        h = radius = 0;
        position = Vec();
        emission = Vec();
        color = Vec();
        refl = Refl::SPEC;
    };
    Object(double rad_, Vec p_, Vec e_, Vec c_, Refl refl_) : radius(rad_), position(p_), emission(e_), color(c_), refl(refl_) {}
    Object(double h_, double rad_, Vec p_, Vec e_, Vec c_, Refl refl_) : h(h_), radius(rad_), position(p_), emission(e_), color(c_), refl(refl_) {}
    double h;      // 高度
    double radius; // 底面半径或者球体半径
    Vec position;  // 中心坐标
    Vec emission;  // 光源强度
    Vec color;     // 颜色
    Refl refl;     // reflection type (DIFFuse, SPECular, REFRactive)
    virtual double intersect(const Ray &r) const = 0;
    virtual Vec compN(const Vec &pos) const = 0; // 计算点pos处的法向量
};
class Sphere : public Object // 场景中的球体
{
public:
    // 球体半径，球心位置，光源，颜色
    //  radius,position,emission,color
    //  refl:reflection type (DIFFuse, SPECular, REFRactive)
    Sphere() : Object() {}
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl refl_) : Object(rad_, p_, e_, c_, refl_) {}
    virtual double intersect(const Ray &r) const
    {
        // returns distance, 0 if nohit Solve
        // 几何关系有|op-td|=r，两边平方后整理出一元二次方程，求解delata，再求解出t的值
        // 方程具体为：t^2 - 2 op·d + op·op - r^2 = 0
        // delta = B^2-4AC = (op·d)^2 - op·op + r^2
        // t = -B(+-)sqrt(delta) / 2A = op·d (+-) sqrt(delta)
        auto op = position - r.o;    // 向量op
        constexpr double eps = 1e-4; // 误差控制
        auto b = op.dot(r.d);        // op·d,也可以表达op和d的余弦值（这里只是记录点积）
        auto det = b * b - op.dot(op) + radius * radius;
        if (det < 0)
            return 0; // 方程无解，返回0表示无交点
        det = sqrt(det);
        auto t = b - det;
        if (t > eps)
            return t;
        // 这里先判断近处的交点，一根射线可能和球有两个交点
        t = b + det;
        if (t > eps)
            return t;
        return 0;
    }
    virtual Vec compN(const Vec &pos) const
    {
        return (pos - position).norm();
    }
};
class Cylinder : public Object // 场景中的圆柱体
{
public:
    Vec topp; // 顶面圆心坐标
    Vec botp; // 底面圆心坐标
    Cylinder() : Object()
    {
        topp = botp = Vec();
    }
    Cylinder(double h_, double rad_, Vec p_, Vec e_, Vec c_, Refl refl_) : Object(h_, rad_, p_, e_, c_, refl)
    {
        topp = Vec(position.x, position.y + h / 2, position.z);
        botp = Vec(position.x, position.y - h / 2, position.z);
    }
    virtual double intersect(const Ray &r) const
    {
        // 圆柱体没有球体那么好的性质，求交点需要分为两步，第一步是判断是否和侧面相交，第二步判断是否和两个底面相交
        // 侧面相交考虑联立射线方程和圆柱侧面方程（无限高度），求解出交点后判断交点是否位于两个底面坐标范围之内
        // 底面考虑和圆的方程联立（x^2+z^2=r^2）
        double m = r.o.x - position.x;
        double n = r.o.z - position.z;
        double A = r.d.x * r.d.x + r.d.z * r.d.z;
        double B = 2 * m * r.d.x + 2 * n * r.d.z;
        double C = m * m + n * n - radius * radius;
        double det = B * B - 4 * A * C;

        if (det >= 0 && A > eps) // 侧面
        {
            det = sqrt(det);
            double t = (-B - det) / (2 * A);
            double yj = r.o.y + t * r.d.y;
            if (t > eps && yj <= topp.y && yj >= botp.y)
                return t;
            // 这里先判断近处的交点，一根射线可能和柱面有两个交点
            t = (-B + det) / (2 * A);
            yj = r.o.y + t * r.d.y;
            if (t > eps && yj <= topp.y && yj >= botp.y)
                return t;
        }

        // 上底面相交判断
        if (r.d.y <= eps)
        {
            double t = (topp.y - r.o.y) / r.d.y;
            double xt = r.o.x + t * r.d.x;
            double zt = r.o.z + t * r.d.z;
            if (t > eps)
                if ((xt - topp.x) * (xt - topp.x) + (zt - topp.z) * (zt - topp.z) <= radius * radius)
                    return t;
        }
        // 下底面判断
        if (r.d.y <= eps){
            double t = (botp.y - r.o.y) / r.d.y;
            double xt = r.o.x + t * r.d.x;
            double zt = r.o.z + t * r.d.z;
            if (t <= eps)
                return 0;
            if ((xt - botp.x) * (xt - botp.x) + (zt - botp.z) * (zt - botp.z) <= radius * radius)
                return t;
        }
        return 0;
    }
    virtual Vec compN(const Vec &pos) const
    {
        if (fabs(pos.y - topp.y) <= eps)
            return Vec(0, -1, 0);
        else if (fabs(pos.y - botp.y) <= eps)
            return Vec(0, 1, 0);
        else
            return Vec(2*pos.x-2*position.x,0, 2*pos.z-2*position.z).norm();
    }
};

inline bool intersect(const Ray &r, Object *&obj, double &t);
Vec radiance(const Ray &r, int depth, unsigned short *Xi);

// 限制x在0~1范围内
inline double clamp(double x) { return x < 0 ? 0 : (x > 1 ? 1 : x); }

inline int toInt(double x) // 伽马矫正，最后生成图像转成RGB值使用
{
    return static_cast<int>(pow(clamp(x), 1 / 2.2) * 255 + .5);
}

inline double erand48(unsigned short *Xi) { return (double)rand() / RAND_MAX; }

inline Vec compNl(const Vec &n, const Ray &r) // 主要区分光线向内还是向外
{
    if (n.dot(r.d) < 0)
        return n;
    return n * -1;
}

inline double computeNnt(bool into, double nc, double nt)
{
    if (into)
        return nc / nt;
    return nt / nc;
}

inline Vec computeTdir(const Ray &r, double nnt, const Vec &n, bool into, double ddn, double cos2t)
{
    return (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
}

Vec radiance(const Ray &r, int depth, unsigned short *Xi)
{
    if (depth > 10)
        return Vec();      // 深度太深直接返回
    double t;              // distance to intersection
    Object *obj = nullptr; // intersected object
    if (!intersect(r, obj, t))
    {
        // std::cout << "nohit\n";
        return {}; // 没有相交，本光线追踪结束，返回黑色
    }
    Vec x = r.o + r.d * t; // 交点
    Vec n = obj->compN(x); // 交点处的法向量
    // Vec n = (x - obj->position).norm();
    Vec nl = compNl(n, r);

    auto maxc = obj->color.greatest();
    Vec f = obj->color;
    if (++depth > 5) // 以一定概率终止递归
    {
        if (erand48(Xi) < maxc)
            f = f * (1 / maxc); // 操作过后颜色值会变大，其中颜色最大的分量会变为1
        else
            return obj->emission; // 递归次数已经较深了，但是像素最大颜色值仍然比较小，直接退出，没有必要再递归下去
        // 如果是光源，则返回光源强度，否则证明是比较暗的场景，返回值也会是0
    }
    // 计算当前颜色贡献
    if (obj->refl == Refl::DIFF) // 漫反射，半球面积分，法线为中心的半球面均匀采样
    // Lambert 模型
    {                                       // Ideal DIFFUSE reflection
        double r1 = 2 * M_PI * erand48(Xi); // 随机一个0~2PI的方向
        double r2 = erand48(Xi);
        double r2s = sqrt(r2);
        // 随机出一组基u,v,nl，在该空间里随机组合产生一个新方向进行递归
        Vec u = ((fabs(nl.x) > 0.1 ? Vec(0, 1, 0) : Vec(1, 0, 0)) % nl).norm();
        // nl的x比较小就用(1,0,0)向量和nl叉乘产生一个新向量u，这样保证(1,0,0)和nl是线性无关的(nl的x分量接近0)
        // 反之可以用(0,1,0)和nl叉乘产生u，主要是保证无关性，(0,1,0)应该只是为了计算方便
        Vec v = nl % u; // 产生一个和u,nl都垂直的v
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + nl * sqrt(1 - r2)).norm();
        return obj->emission + f.mult(radiance({x, d}, depth, Xi));
    }
    else if (obj->refl == Refl::SPEC) // 镜面反射
    {                                 // Ideal SPECULAR reflection
        return obj->emission + f.mult(radiance({x, r.d - n * 2 * n.dot(r.d)}, depth, Xi));
    }

    // 折射处理
    Ray reflRay = {x, r.d - n * 2 * n.dot(r.d)};    // 反射射线
    bool into = n.dot(nl) > 0;                      // 内部到外部还是外部到内部
    constexpr double nc = 1, nt = 1.5;              // 空气、玻璃折射率
    double nnt = computeNnt(into, nc, nt);          // 折射率比值，折射前/折射后，外部光线是nc/nt，反之相反
    double ddn = r.d.dot(nl);                       // 法线、光线余弦值
    double cos2t = 1 - nnt * nnt * (1 - ddn * ddn); // 主要是折射定律推导，判断出射角是否大于90度，cost^2肯定是大于0的
    // 但是当出现全反射时，折射定律是不满足的，因此这样计算的cost^2就会变得小于0
    if (cos2t < 0)
        return obj->emission + f.mult(radiance(reflRay, depth, Xi)); // 内部射出光线全反射的情况
    // 正常处理折射、反射的情况
    Vec tdir = computeTdir(r, nnt, n, into, ddn, cos2t); // 处理折射方向
    constexpr double a = nt - nc;
    constexpr double b = nt + nc;
    constexpr double R0 = a * a / (b * b);
    // 菲涅尔项(Fresnel Term): 反射率取决于入射角度，入射光与法线的夹角越大
    // 反射的能量越多，折射越少，这里是近似计算
    // R(theta) = R0 + (1-R0)*(1-cos(theta))^5,R0=((n1-n2)/(n1+n2))^2,n1,n2代表折射率
    double c = 1 - (into ? -ddn : tdir.dot(n));
    double Re = R0 + (1 - R0) * c * c * c * c * c; // 反射部分
    double Tr = 1 - Re;                            // 折射部分
    double P = .25 + .5 * Re;
    double RP = Re / P;
    double TP = Tr / (1 - P);
    // 深度小于3的时候分别计算折射、反射两条光线，深度大的时候，只计算一条光线，并根据光线能量乘占比再除以概率恢复
    // 两条光线的效果。在这里以一定概率选择计算反射部分或者折射部分
    return obj->emission + f.mult(depth > 2 ? (erand48(Xi) < P ? radiance(reflRay, depth, Xi) * RP : radiance({x, tdir}, depth, Xi) * TP)
                                            : radiance(reflRay, depth, Xi) * Re + radiance({x, tdir}, depth, Xi) * Tr);
}

std::vector<Object *> objs; // 存储场景

inline bool intersect(const Ray &r, Object *&obj, double &t)
{
    constexpr double inf = 1e20;
    t = inf;
    for (auto &it : objs)
    {
        auto d = it->intersect(r);
        if (d && d < t)
        {
            t = d;
            obj = it;
        }
    }
    return t < inf;
}

void initialization()
{
    objs.push_back(new Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), Refl::DIFF));   // Left
    objs.push_back(new Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), Refl::DIFF)); // Rght
    objs.push_back(new Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), Refl::DIFF));         // Back
    objs.push_back(new Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), Refl::DIFF));               // Frnt
    objs.push_back(new Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), Refl::DIFF));         // Botm
    objs.push_back(new Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), Refl::DIFF)); // Top

    objs.push_back(new Sphere(8, Vec(27, 8, 47), Vec(), Vec(1, 1, 1) * .999, Refl::SPEC));           // Mirr
    objs.push_back(new Sphere(8, Vec(73, 8, 78), Vec(), Vec(1, 1, 1) * .999, Refl::REFR));           // Glass
    objs.push_back(new Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(12, 12, 12), Vec(), Refl::DIFF)); // Lite
    objs.push_back(new Cylinder(10, 4, Vec(36, 5, 86), Vec(), Vec(1, 1, 1) * .999, Refl::DIFF));     // Cylinder

    objs.push_back(new Cylinder(20, 1, Vec(6, 8, 76), Vec(), Vec(.4, .6, .3) * .999, Refl::DIFF));          // Cylinder
    objs.push_back(new Sphere(5, Vec(6, 23, 76), Vec(5, 5, 5), Vec(1.0, 0.842, 0.706) * .999, Refl::DIFF)); // lite2
}

void initialization2()
{
    Vec Cen(50, 40.8, -860);
    objs.push_back(new Sphere(1600, Vec(1, 0, 2) * 3000, Vec(1, .9, .8) * 1.2e1 * 1.56 * 2, Vec(), Refl::DIFF));  // sun
    objs.push_back(new Sphere(1560, Vec(1, 0, 2) * 3500, Vec(1, .5, .05) * 4.8e1 * 1.56 * 2, Vec(), Refl::DIFF)); // horizon sun2
    // objs.push_back(new Sphere(10000,Cen+Vec(0,0,-200), Vec(0.0627, 0.188, 0.569)*6e-2*8, Vec(.7,.7,1)*.25,Refl::DIFF));// sky
    objs.push_back(new Sphere(10000, Cen + Vec(0, 0, -200), Vec(0.00063842, 0.02001478, 0.28923243) * 6e-2 * 8, Vec(.7, .7, 1) * .25, Refl::DIFF)); // sky
    objs.push_back(new Sphere(100000, Vec(50, -100000, 0), Vec(), Vec(.3, .3, .3), Refl::DIFF));                                                    // grnd
    objs.push_back(new Sphere(110000, Vec(50, -110048.5, 0), Vec(.9, .5, .05) * 4, Vec(), Refl::DIFF));                                             // horizon brightener
    // objs.push_back(new Sphere(4e4, Vec(50, -4e4-30, -3000),  Vec(),Vec(.2,.2,.2),Refl::DIFF));// mountains
    objs.push_back(new Sphere(3.99e4, Vec(50, -3.99e4 + 20.045, -3000), Vec(), Vec(.7, .7, .7), Refl::DIFF)); // mountains snow
    objs.push_back(new Sphere(26.5, Vec(22, 26.5, 42), Vec(), Vec(1, 1, 1) * .596, Refl::SPEC));              // white Mirr
    objs.push_back(new Sphere(13, Vec(75, 13, 82), Vec(), Vec(.96, .96, .96) * .96, Refl::REFR));             // Glas
    objs.push_back(new Sphere(22, Vec(87, 22, 24), Vec(), Vec(.6, .6, .6) * .696, Refl::REFR));  // Glas2
    objs.push_back(new Cylinder(10, 4, Vec(36, 5, 86), Vec(), Vec(.5, .5, .5) * .999, Refl::DIFF));  // Cylinder
}

int main(int argc, char *argv[])
{
    constexpr int width = 1024;
    constexpr int height = 768;
    int samples = argc == 2 ? atoi(argv[1]) / 4 : 1; // # samples
    // 一个像素采样四个格，这样可以保证采样在该像素分布尽量分散，不是集中在某一点
    Ray camera{{50, 52, 295.6}, Vec{0, -0.042612, -1}.norm()}; // camera pos, dir
    initialization();
    // initialization2();
    Vec cx = Vec{width * .5135 / height};    // cx,cy完成对成像平面的量化，即遍历时的步长
    Vec cy = (cx % camera.d).norm() * .5135; // 0.5135实际上是给出FOV（视场角）
    Vec r;                                   // 每次采样的颜色
    std::vector<Vec> c(width * height);      // 图像

#pragma omp parallel for schedule(dynamic, 1) private(r) // OpenMP

    for (int y = 0; y < height; y++)
    { // Loop over image rows
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samples * 4, 100. * y / (height - 1));
        // Loop over cols
        for (unsigned short x = 0, Xi[3] = {0, 0, static_cast<unsigned short>(y * y * y)}; x < width; x++)
            for (int sy = 0, i = (height - y - 1) * width + x; sy < 2; sy++) // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec())
                { // 2x2 subpixel cols
                    for (int s = 0; s < samples; s++)
                    {
                        double r1 = 2 * erand48(Xi);
                        double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        // 让采样像素中间的概率更大，类似一种滤波操作
                        double r2 = 2 * erand48(Xi);
                        double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / width - .5) +
                                cy * (((sy + .5 + dy) / 2 + y) / height - .5) + camera.d;
                        auto mult140 = camera.o + d * 140;                             // 限定相机原点与成像平面距离140个单位
                        r = r + radiance({mult140, d.norm()}, 0, Xi) * (1. / samples); // 蒙特卡洛积分
                    }                                                                  // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + Vec{clamp(r.x), clamp(r.y), clamp(r.z)} * .25;
                    // 采样四个格，最后取颜色平均值，直接计算时乘1/4
                }
    }
    FILE *f = fopen("image.ppm", "w"); // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
    for (const auto &vec : c)
        fprintf(f, "%d %d %d ", toInt(vec.x), toInt(vec.y), toInt(vec.z));
}
