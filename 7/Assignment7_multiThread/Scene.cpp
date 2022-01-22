//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"

void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// // Implementation of Path Tracing
// // Randomly choose ONE dirction wi~pdf(w)
// Vector3f Scene::castRay(const Ray &ray, int depth) const
// {
//     // TO DO Implement Path Tracing Algorithm here
//     if (depth > this->maxDepth)
//     {
//         return Vector3f(0.0, 0.0, 0.0);
//     }

//     Vector3f L_dir(0.0, 0.0, 0.0);
    

//     //Shot a ray from eye_pos to p(eye_pos to reflect)
//     // p, need to intersect
//     Intersection intersection1 = Scene::intersect(ray);
//     if (!intersection1.happened)
//     {
//         return Vector3f(0.0, 0.0, 0.0);
//     }
//     if (intersection1.m->hasEmission())
//     {
//         return intersection1.m->getEmission();
//     }
    
//     Vector3f p = intersection1.coords;

//     // shot a ray form p to x'(reflect to light)
//     // Sample the light
//     // Random sample one light at x'
//     // pass reference parameters
//     Intersection intersection;
//     float pdf;
//     sampleLight(intersection, pdf);  // intersection at x'(xx), pdf = 1 / A
//     Vector3f xx = intersection.coords;  // x'
//     Vector3f dir = p - xx;
//     dir = normalize(dir);
//     Ray ray1(p, dir);
//     Intersection intersection2 = Scene::intersect(ray1);
//     // If the ray is not blocked in the middle
//     // if (intersection2.m == nullptr)
//     // {
//     //     std::cout << "xxxxxxxxxx" << std::endl;
//     // }
//     // if (intersection2.m->getType() != DIFFUSE)
//     // {
//     //     // L_dir = Li * fr * cos(theta) * cos(theta') / |x' - p|^2 / pdf_light
//     //     Vector3f Li = intersection2.m->getEmission();
//     //     // eval wi wo N
//     //     Vector3f fr1 = intersection2.m->eval(ray.direction, ray1.direction, intersection1.normal);
//     //     float costheta = dotProduct(intersection1.normal.normalized(), ray1.direction);
//     //     float costhetax = dotProduct(intersection2.normal.normalized(), ray1.direction);
//     //     float distance_xx_p_2 =  std::pow(xx.x - p.x, 2) + std::pow(xx.y - p.y, 2) + std::pow(xx.z - p.z, 2);
//     //     L_dir += (Li * fr1 * costheta * costhetax / distance_xx_p_2 / pdf); 
//     // }
//     if (intersection2.distance - dir.norm() > -EPSILON)
//     {
//         // L_dir = Li * fr * cos(theta) * cos(theta') / |x' - p|^2 / pdf_light
//         Vector3f Li = intersection2.m->getEmission();
//         // eval wi wo N
//         Vector3f fr1 = intersection2.m->eval(ray.direction, ray1.direction, intersection1.normal);
//         float costheta = dotProduct(intersection1.normal.normalized(), ray1.direction);
//         float costhetax = dotProduct(intersection2.normal.normalized(), ray1.direction);
//         float distance_xx_p_2 =  std::pow(xx.x - p.x, 2) + std::pow(xx.y - p.y, 2) + std::pow(xx.z - p.z, 2);
//         L_dir += (Li * fr1 * costheta * costhetax / distance_xx_p_2 / pdf); 
//     }

//     Vector3f L_indir(0.0, 0.0, 0.0);
//     // Russian Roulette
//     srand(time(NULL));
//     float ksi = rand() % (1000) / float(1000);
//     if (ksi > RussianRoulette)
//     {
//         // nothing to do
//         // L_indir = 0
//     }
//     else
//     {
//         // Uniform sample the hemisphere toward wi (pdf_hemi = 1 / 2pi)
//         Vector3f wo = intersection1.m->sample(ray.direction, intersection1.normal);
//         // Trace a ray r(p, wo)
//         Ray ray2(p, wo);
//         // from shader point to another reflect point(next hitpoint, not light)
//         Intersection intersection3 = Scene::intersect(ray2);
//         // If the ray hit a non-emitting object at q
//         if (intersection3.m->getType() == DIFFUSE)
//         {
//             // L_indir = shade(q, -wi) * fr * cos(theta) / pdf_hemi / P_RR
//             // theta
//             // intersection3.normal
//             // another dir
//             // fr is linked to shader point
//             Vector3f fr2 = intersection1.m->eval(ray.direction, ray2.direction, intersection1.normal);
//             float costheta = dotProduct(intersection1.normal.normalized(), ray2.direction);
//             Ray ray3(intersection3.coords, ray2.direction.normalized());
//             L_indir += (castRay(ray3, depth + 1) * fr2 * costheta / (1 / 2 * M_PI) / RussianRoulette);
//         }
//     }

//     return L_dir + L_indir;
// }

// Implementation of Path Tracing
// Randomly choose ONE dirction wi~pdf(w)
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    Vector3f L_dir(0.0, 0.0, 0.0);
    Vector3f L_indir(0.0, 0.0, 0.0);

    Intersection intersection = Scene::intersect(ray);
    // exceed the scene
    if (!intersection.happened)
    {
        return Vector3f(0.0, 0.0, 0.0);
    }
    // see toward to light directly
    if (intersection.m->hasEmission())
    {
        if (depth == 0)
            return intersection.m->m_emission;
        else 
            return Vector3f(0.0, 0.0, 0.0);
    }

    // L_dir
    // sample the light
    float light_pdf;
    Intersection light_inter;
    sampleLight(light_inter, light_pdf);
    // whether the ray is blocked in the middle 
    Vector3f p = intersection.coords;
    Vector3f xx = light_inter.coords;  // x'
    // Vector3f dir = p - xx;  // error
    Vector3f dir = xx - p;
    dir = normalize(dir);
    Ray ray1(p, dir);
    Intersection whetherLightBlocked = Scene::intersect(ray1);
    // distance
    // float distance1_xx_p_pow = dir.x * dir.x + dir.y *dir.y + dir.z * dir.z;  // error dir has been normalized
    float distance_xx_p_pow = (xx - p).x * (xx - p).x + (xx - p).y * (xx - p).y +
                              (xx - p).z * (xx - p).z;

    const float epsilon = 0.0005f;
    // if (std::fabs(distance1_xx_p - distance2_p_whetherLightBlocked) < EPSILON)  // Black screen
    if ((xx - whetherLightBlocked.coords).norm() < 1e-2)
    {
        // L_dir = Li * fr * cos(theta) * cos(theta') / |x' - p|^2 / pdf_light
        Vector3f Li = light_inter.emit;
        Vector3f fr = intersection.m->eval(ray.direction, ray1.direction, intersection.normal);
        float costheta = dotProduct(intersection.normal, ray1.direction);
        // float costhetax = dotProduct(light_inter.normal, ray1.direction);  // error
        float costhetax = dotProduct(light_inter.normal, -ray1.direction);
        L_dir = Li * fr * costheta * costhetax / distance_xx_p_pow / light_pdf;
    }

    // L_indir
    // Russian Roulette
    // srand(time(NULL));
    // float ksi = rand() % (1000) / float(1000);
    if (get_random_float() > RussianRoulette)
    {
        // L_indir nothing to do
        return L_dir;
    }
    else
    {
        // Uniform sample the hemisphere toward wi (pdf_hemi = 1 / 2 * pi)
        Vector3f l_exit_dir = intersection.m->sample(ray.direction, intersection.normal);
        l_exit_dir = normalize(l_exit_dir);
        Ray ray2(p, l_exit_dir);
        Intersection inter_hemisphere = Scene::intersect(ray2);

        // If the ray hit a non-emitting object at q
        if (inter_hemisphere.happened && !inter_hemisphere.m->hasEmission())
        {
            Vector3f q = inter_hemisphere.coords;
            // L_indir = shade(q, -wi) * fr * cos(theta) / pdf_hemi / P_RR
            // float pdf_error = inter_hemisphere.m->pdf(ray.direction, l_exit_dir, intersection.normal);  // error
            float pdf = intersection.m->pdf(ray.direction, ray2.direction, intersection.normal);
            Vector3f fr = intersection.m->eval(ray.direction, ray2.direction, intersection.normal);
            float costheta = dotProduct(intersection.normal, ray2.direction);
            L_indir = castRay(ray2, depth + 1) * fr * costheta / pdf / RussianRoulette;
        }
    }
    
    return L_dir + L_indir;
	// Intersection inter = intersect(ray);

	// if (inter.happened)
	// {
	// 	// 如果射线第一次打到光源，直接返回
	// 	if (inter.m->hasEmission())
	// 	{
	// 		if (depth == 0) 
	// 		{
	// 			return inter.m->getEmission();
	// 		}
	// 		else return Vector3f(0, 0, 0);
	// 	}

	// 	Vector3f L_dir(0, 0, 0);
	// 	Vector3f L_indir(0, 0, 0);

	// 	// 随机 sample 灯光，用该 sample 的结果判断射线是否击中光源
	// 	Intersection lightInter;
	// 	float pdf_light = 0.0f;
	// 	sampleLight(lightInter, pdf_light);

	// 	// 物体表面法线
	// 	auto& N = inter.normal;
	// 	// 灯光表面法线
	// 	auto& NN = lightInter.normal;

	// 	auto& objPos = inter.coords;
	// 	auto& lightPos = lightInter.coords;

	// 	auto diff = lightPos - objPos;
	// 	auto lightDir = diff.normalized();
	// 	float lightDistance = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;

	// 	Ray light(objPos, lightDir);
	// 	Intersection light2obj = intersect(light);

	// 	// 如果反射击中光源
	// 	if (light2obj.happened && (light2obj.coords - lightPos).norm() < 1e-2)
	// 	{
	// 		Vector3f f_r = inter.m->eval(ray.direction, lightDir, N);
	// 		L_dir = lightInter.emit * f_r * dotProduct(lightDir, N) * dotProduct(-lightDir, NN) / lightDistance / pdf_light;
	// 	}

	// 	if (get_random_float() < RussianRoulette)
	// 	{
	// 		Vector3f nextDir = inter.m->sample(ray.direction, N).normalized();

	// 		Ray nextRay(objPos, nextDir);
	// 		Intersection nextInter = intersect(nextRay);
	// 		if (nextInter.happened && !nextInter.m->hasEmission())
	// 		{
	// 			float pdf = inter.m->pdf(ray.direction, nextDir, N);
	// 			Vector3f f_r = inter.m->eval(ray.direction, nextDir, N);
	// 			L_indir = castRay(nextRay, depth + 1) * f_r * dotProduct(nextDir, N) / pdf / RussianRoulette;
	// 		}
	// 	}

	// 	return L_dir + L_indir;
	// }

	// return Vector3f(0, 0, 0);
}