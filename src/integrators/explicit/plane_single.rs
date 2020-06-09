use std::f32::consts::PI;
use std::ops::Add;
use std::process::exit;

use crate::accel::*;
use crate::geometry::Mesh;
use crate::integrators::*;
use crate::math::*;
use crate::samplers;
use crate::structure::AABB;
use crate::volume::*;
use cgmath::{ElementWise, EuclideanSpace, InnerSpace, Point2, Point3, Vector3};
use nalgebra::ComplexField;
use rand::Rng;
use rand::thread_rng;
use rand::distributions::{Distribution, WeightedIndex};
use rand_distr::Normal as EasyNormal;
use statrs::distribution::ContinuousCDF;
use statrs::distribution::Normal;

fn clamp<T: PartialOrd>(v: T, min: T, max: T) -> T {
    if v < min {
        min
    } else if v > max {
        max
    } else {
        v
    }
}

fn direct_int(a: f32, b: f32, c: f32,A: f32,B: f32) -> f32 {
    if c<f32::EPSILON
    {
        return 0.0;
    };
    let t = f32::atan2(B, -A);
    let t1 = f32::atan2(c, b);
    let t2 = f32::atan2(c, -a);
    let sign1 = if t < t1 { -1.0 } else { 1.0 };
    let sign2 = if t1 <= t && t <= t2 { -1.0 } else { sign1 };
    let int1 = {
        if b>f32::EPSILON{
            if t < t1 {
                b * (A * f32::ln(f32::cos(t1) / f32::cos(t).powi(2)) + B * (2.0 * t - t1))
            } else {
                b * (A * f32::ln(1.0 / f32::cos(t1)) + B * t1)
            }
        }
        else{
            0.0
        }
    };
    //assert!(int1>=0.0);
    let int2 ={
        if t1 <= t && t <= t2 {
            c * (A * (2.0 * t - t1 - t2) + B * f32::ln(f32::sin(t).powi(2) / (f32::sin(t1) * f32::sin(t2))))
        } else {
            c * (A * (t2 - t1) + B * f32::ln(f32::sin(t2) / f32::sin(t1))) * sign1
        }
    };
    //assert!(int2>=0.0);
    let int3 = {
        if a>f32::EPSILON{
            if t > t2 {
                //sign = -1.0;
                -a * (A * f32::ln(-f32::cos(t2) / f32::cos(t).powi(2)) + B * (2.0 * t - t2 - PI))
                
            } else {
                -a * (A * f32::ln(-f32::cos(t2)) + B * (PI - t2)) * sign2
            }  
        }
        else {
            0.0
        }
    };
    //assert!(int3>=0.0);
    (int1 + int2 + int3) / PI
}

fn generate_sin_sample(bias: f32) -> f32 {
    let mut rng = rand::thread_rng();
    let p: f32 = rng.gen(); // 生成 [0, 1] 区间内的随机数
    let sample = (1.0 / PI) * f32::acos(1.0 - 2.0 * p) - bias / PI;
    sample.rem_euclid(1.0) // 确保结果在 [0, 1] 区间内
}

fn sin_pdf(x: f32,bias:f32) -> f32 {
    (PI/2.0*f32::sin(PI*x+bias)).abs()
}

fn generate_normal_sample(std_dev: f32) -> f32 {
    let normal = EasyNormal::new(0.5, std_dev as f32).unwrap();
    let mut rng = thread_rng();

    loop {
        let sample = rng.sample(normal) as f32;
        if (0.0..=1.0).contains(&sample) {
            return sample;
        }
    }
}

fn normal_pdf(x: f32,std_dev:f32,normalization_factor:f32) -> f32 {
    let mean = 0.5;
    let exponent = -(x - mean).powi(2) / (2.0 * std_dev.powi(2));
    (1.0 / (std_dev * (2.0 * std::f32::consts::PI).sqrt())) * exponent.exp() / normalization_factor
}

// 采样点与角度转换
fn calculate_angle_with_ab_corrected_v2(x_i: f32, y_i: f32, a: f32, b: f32) -> f32 {
    let dx = x_i - a;
    let dy = y_i - b;
    let angle = dy.atan2(dx) / std::f32::consts::PI;
    angle.rem_euclid(1.0)
}

// 计算函数 f(m, n)
fn calculate_f(m: f32, n: f32) -> f32 {
    if n == 0.0 || m == 0.0 {
        return 0.0;
    }
    m * ((1.0 + (n / m).powi(2)).sqrt() + n / m).ln()
}

fn generate_len_sample(u: f32, v: f32, a: f32, b: f32) -> f32 {
    // 四个小矩形的长宽
    let rectangles = [
        (a, b), (u - a, b), (a, v - b), (u - a, v - b), 
        (b, a), (b, u - a), (v - b, a), (v - b, u - a)
    ];

    // 计算每个小矩形的 f(m, n)
    let mut f_values = Vec::new();
    for (i, &(m, n)) in rectangles.iter().enumerate() {
        f_values.push((calculate_f(m, n), i));
    }

    let total_f: f32 = f_values.iter().map(|f| f.0).sum();
    let weights: Vec<_> = f_values.iter().map(|f| f.0 / total_f).collect();
    
    let dist = match WeightedIndex::new(&weights) {
        Ok(d) => d,
        Err(e) => {
            // 打印错误信息和权重数组
            eprintln!("Weights: {:?}", weights);
            eprintln!("f_values: {:?}", f_values);
            eprintln!("rectangles: {:?}", rectangles);
            exit(0);
        }
    };
    
    let mut rng = thread_rng();
    let sampled_index = dist.sample(&mut rng);
    let (c_a, c_b) = rectangles[sampled_index];
    let y: f32 = rand::random(); // 生成一个随机数

    let t = ((c_a.powi(2) + c_b.powi(2)).sqrt() + c_b) / c_a;
    let n = t.powf(y);
    let x = c_a * (n - 1.0 / n) / 2.0;

    // 根据索引确定采样点
    let (sampled_x, sampled_y) = match sampled_index {
        0..=3 => {
            let sampled_x = if sampled_index % 2 == 0 { 0.0 } else { u };
            let sampled_y = if sampled_index < 2 { b - x } else { b + x };
            (sampled_x, sampled_y)
        },
        _ => {
            let sampled_x = if sampled_index % 2 == 0 { a - x } else { a + x };
            let sampled_y = if sampled_index < 6 { 0.0 } else { v };
            (sampled_x, sampled_y)
        },
    };
    calculate_angle_with_ab_corrected_v2(sampled_x, sampled_y, a, b)
}

fn len_pdf(intergrate:f32,length:f32) ->f32 {
    length/intergrate
}

#[derive(PartialEq, Clone, Debug)]
pub enum PlaneType {
    UV,
    VT,
    UT,
    UAlphaT,
}

// Helper on the light source
pub struct RectangularLightSource {
    pub o: Point3<f32>,
    pub n: Vector3<f32>,
    pub u: Vector3<f32>,
    pub v: Vector3<f32>,
    pub u_l: f32,
    pub v_l: f32,
    pub emission: Color,
}
impl RectangularLightSource {
    pub fn from_shape(emitter: &Mesh) -> Self {
        info!("Emitter vertices: {:?}", emitter.vertices);
        info!("Emitter indices: {:?}", emitter.indices);
        if emitter.vertices.len() != 3 && emitter.indices.len() != 2 {
            panic!("Only support rectangular emitters");
        }
        let o = Point3::from_vec(emitter.vertices[0]);
        let u = emitter.vertices[1] - emitter.vertices[0];
        let v = emitter.vertices[3] - emitter.vertices[0];
        let u_l = u.magnitude();
        let v_l = v.magnitude();

        // Normalize vectors
        let u = u / u_l;
        let v = v / v_l;

        info!("o : {:?}", o);
        info!("u : {:?} ({})", u, u_l);
        info!("v : {:?} ({})", v, v_l);

        let n = u.cross(v);
        info!("Light source normal: {:?}", n);

        let emission = match emitter.emission {
            crate::geometry::EmissionType::Color { v } => v,
            _ => panic!(),
        };

        RectangularLightSource {
            o,
            n,
            u,
            v,
            u_l,
            v_l,
            emission,
        }
    }
}

// TODO: The class is very similar to photon planes
#[derive(Clone, Debug)]
pub struct SinglePhotonPlane {
    o: Point3<f32>,            // Plane origin
    pub d0: Vector3<f32>,      // First edge (normalized)
    pub d1: Vector3<f32>,      // Second edge (normalized)
    length0: f32,              // First edge length
    length1: f32,              // Second edge length
    pub sample: Point2<f32>,   // Random number used to generate the plane (TODO: Unused?)
    pub plane_type: PlaneType, // How the plane have been generated
    pub weight: Color,         // This factor will vary between the different light sources
    id_emitter: usize,
    pub sample_alpha: f32,
}

#[derive(Debug)]
pub struct PhotonPlaneIts {
    pub t_cam: f32,
    t0: f32,
    t1: f32,
    inv_det: f32,
}
impl BVHElement<PhotonPlaneIts> for SinglePhotonPlane {
    fn aabb(&self) -> AABB {
        let p0 = self.o + self.d0 * self.length0;
        let p1 = self.o + self.d1 * self.length1;
        let p2 = p0 + self.d1 * self.length1;
        let mut aabb = AABB::default();
        aabb = aabb.union_vec(&self.o.to_vec());
        aabb = aabb.union_vec(&p0.to_vec());
        aabb = aabb.union_vec(&p1.to_vec());
        aabb = aabb.union_vec(&p2.to_vec());
        aabb
    }
    // Used to construct AABB (by sorting elements)
    fn position(&self) -> Point3<f32> {
        // Middle of the photon plane
        // Note that it might be not ideal....
        self.o + self.d0 * self.length0 * 0.5 + self.d1 * self.length1 * 0.5
    }
    // This code is very similar to triangle intersection
    // except that we loose one test to make posible to
    // intersect planar primitives
    fn intersection(&self, r: &Ray) -> Option<PhotonPlaneIts> {
        let e0 = self.d0 * self.length0;
        let e1 = self.d1 * self.length1;

        let p = r.d.cross(e1);
        let det = e0.dot(p);
        if det.abs() < 1e-5 {
            return None;
        }

        let inv_det = 1.0 / det;
        let t = r.o - self.o;
        let t0 = t.dot(p) * inv_det;
        if t0 < 0.0 || t0 > 1.0 {
            return None;
        }

        let q = t.cross(e0);
        let t1 = r.d.dot(q) * inv_det;
        if t1 < 0.0 || t1 > 1.0 {
            return None;
        }

        let t_cam = e1.dot(q) * inv_det;
        if t_cam <= r.tnear || t_cam >= r.tfar {
            return None;
        }

        // Scale to the correct distance
        // In order to use correctly transmittance sampling
        let t1 = t1 * self.length1;
        let t0 = t0 * self.length0;

        Some(PhotonPlaneIts {
            t_cam,
            t0,
            t1,
            inv_det,
        })
    }
}
impl SinglePhotonPlane {
    pub fn light_position(
        &self,
        light: &RectangularLightSource,
        plane_its: &PhotonPlaneIts,
    ) -> Point3<f32> {
        match self.plane_type {
            PlaneType::UV => light.o + light.u * plane_its.t0 + light.v * plane_its.t1,
            PlaneType::VT | PlaneType::UT | PlaneType::UAlphaT => self.o + self.d0 * plane_its.t0,
        }
    }
    pub fn contrib(&self, d: &Vector3<f32>) -> Color {
        let jacobian = self.d1.cross(self.d0).dot(*d).abs();
        self.weight / jacobian
    }
    pub fn new(
        t: PlaneType,
        light: &RectangularLightSource,
        d: Vector3<f32>,
        sample: Point2<f32>,
        sample_alpha: f32,
        t_sampled: f32,
        id_emitter: usize,
        sigma_s: Color,
    ) -> Self {
        match t {
            PlaneType::UV => {
                let o_plane = light.o + d * t_sampled;
                SinglePhotonPlane {
                    o: o_plane,
                    d0: light.u,
                    d1: light.v,
                    length0: light.u_l,
                    length1: light.v_l,
                    sample,
                    plane_type: PlaneType::UV,
                    // Sigma_s is used to cancel out (as we sample a distance in planes)
                    weight: std::f32::consts::PI * light.emission / sigma_s,
                    id_emitter,
                    sample_alpha,
                }
            }
            PlaneType::VT => {
                let o_plane = light.o + light.u * light.u_l * sample.x;
                SinglePhotonPlane {
                    o: o_plane,
                    d0: light.v,
                    d1: d,
                    length0: light.v_l,
                    length1: t_sampled,
                    sample,
                    plane_type: PlaneType::VT,
                    // f() / (p(w_l) * p(u))
                    weight: std::f32::consts::PI * light.u_l * light.emission,
                    id_emitter,
                    sample_alpha,
                }
            }
            PlaneType::UT => {
                let o_plane = light.o + light.v * light.v_l * sample.y;
                SinglePhotonPlane {
                    o: o_plane,
                    d0: light.u,
                    d1: d,
                    length0: light.u_l,
                    length1: t_sampled,
                    sample,
                    plane_type: PlaneType::UT,
                    // f() / (p(w_l) * p(u))
                    weight: std::f32::consts::PI * light.v_l * light.emission,
                    id_emitter,
                    sample_alpha,
                }
            }
            PlaneType::UAlphaT => {
                // Function to compute intersection on 2D AABB
                let plane2d_its = |d: Vector2<f32>, o: Point2<f32>| {
                    let t_0 = (-o.to_vec()).div_element_wise(d);
                    let t_1 = (Vector2::new(light.u_l, light.v_l) - o.to_vec()).div_element_wise(d);
                    let t_max_coord = Vector2::new(t_0.x.max(t_1.x), t_0.y.max(t_1.y));
                    o + d * t_max_coord.x.min(t_max_coord.y)
                };

                // 1) Generate point on the light
                let o_plane = Point2::new(sample.x * light.u_l, sample.y * light.v_l);

                // 2) Generate orientation of the plane
                let alpha = std::f32::consts::PI * sample_alpha;
                let d_plane: Vector2<f32> = Vector2::new(alpha.cos(), alpha.sin());

                // 3) Compute the intersection point from the 2D light source
                //    to compute the segment on it
                let p1_2d = plane2d_its(d_plane, o_plane);
                let p2_2d = plane2d_its(-d_plane, o_plane);

                // 4) Convert the 2D end-point of the segment in 3D
                let p1 = light.o + p1_2d.x * light.u + p1_2d.y * light.v;
                let p2 = light.o + p2_2d.x * light.u + p2_2d.y * light.v;

                // 5) Compute the segment in 3D
                let u_plane = p2 - p1;
                let u_plane_length = u_plane.magnitude();
                let u_plane = u_plane / u_plane_length;

                SinglePhotonPlane {
                    // Plane geometry
                    o: p1,
                    d0: u_plane,
                    d1: d,
                    // d0 length = segment length
                    length0: u_plane_length,
                    // d1 length = use transmittance sampling
                    length1: t_sampled,
                    // storing the random sample in case
                    sample,
                    // Mark that it is a (u, alpha)-plane
                    plane_type: PlaneType::UAlphaT,
                    // Compute the plane weight
                    // f() / (p(w_l) * p(e)) where
                    //  p(w_l) = cos(..) / PI, the cos cancel out with emission
                    //  p(e) = |e| / A
                    weight: std::f32::consts::PI * light.emission * (light.u_l * light.v_l)
                        / u_plane_length,
                    // Store the emitter ID
                    id_emitter,
                    sample_alpha,
                }
            }
        }
    }
}

#[derive(PartialEq)]
pub enum SinglePlaneStrategy {
    UV,
    VT,
    UT,
    Average,
    DiscreteMIS,
    UAlpha,
    ContinousMIS,
    SMISAll(usize),
    SMISJacobian(usize),
    ProxySample,
}

pub struct IntegratorSinglePlane {
    pub nb_primitive: usize,
    pub strategy: SinglePlaneStrategy,
    pub stratified: bool,
}

impl Integrator for IntegratorSinglePlane {
    fn compute(
        &mut self,
        sampler: &mut dyn Sampler,
        accel: &dyn Acceleration,
        scene: &Scene,
    ) -> BufferCollection {
        if scene.volume.is_none() {
            panic!("Volume integrator need a volume (add -m )");
        }
        // Extract the light source
        let emitters = scene
            .meshes
            .iter()
            .filter(|m| m.is_light())
            .collect::<Vec<_>>();
        let rect_lights = {
            emitters
                .iter()
                .map(|emitter| {
                    info!("Emitter vertices: {:?}", emitter.vertices);
                    info!("Emitter indices: {:?}", emitter.indices);
                    if emitter.vertices.len() != 3 && emitter.indices.len() != 2 {
                        panic!("Only support rectangular emitters");
                    }
                    RectangularLightSource::from_shape(emitter)
                })
                .collect::<Vec<_>>()
        };

        let generate_plane = |t: PlaneType,
                              light: &[RectangularLightSource],
                              id_emitter: usize,
                              sampler: &mut dyn Sampler,
                              m: &HomogenousVolume|
         -> SinglePhotonPlane {
            let d = {
                let mut d_out = cosine_sample_hemisphere(sampler.next2d());
                while d_out.z == 0.0 {
                    // Start to generate a new plane again
                    d_out = cosine_sample_hemisphere(sampler.next2d());
                }
                let frame = Frame::new(light[id_emitter].n);
                frame.to_world(d_out)
            };

            // FIXME: Faking position as it is not important for sampling the transmittance
            let ray_med = Ray::new(light[id_emitter].o, d);
            // TODO: Check if it is the code
            // Need to check the intersection distance iff need to failed ...
            // ray_med.tfar = intersection_distance;
            let mrec = m.sample(&ray_med, sampler.next());

            // Sample planes
            let sample = sampler.next2d();
            SinglePhotonPlane::new(
                t,
                &light[id_emitter],
                d,
                sample,
                sampler.next(),
                mrec.continued_t,
                id_emitter,
                m.sigma_s,
            )
        };

        // Create the planes
        let m = scene.volume.as_ref().unwrap();
        let mut planes = vec![];
        let mut number_plane_gen = 0;
        while planes.len() < self.nb_primitive {
            let id_emitter = (sampler.next() * rect_lights.len() as f32) as usize;
            match self.strategy {
                SinglePlaneStrategy::UT => planes.push(generate_plane(
                    PlaneType::UT,
                    &rect_lights,
                    id_emitter,
                    sampler,
                    m,
                )),
                SinglePlaneStrategy::VT => planes.push(generate_plane(
                    PlaneType::VT,
                    &rect_lights,
                    id_emitter,
                    sampler,
                    m,
                )),
                SinglePlaneStrategy::UV => planes.push(generate_plane(
                    PlaneType::UV,
                    &rect_lights,
                    id_emitter,
                    sampler,
                    m,
                )),
                SinglePlaneStrategy::DiscreteMIS | SinglePlaneStrategy::Average => {
                    // Generate 3 planes
                    planes.push(generate_plane(
                        PlaneType::UV,
                        &rect_lights,
                        id_emitter,
                        sampler,
                        m,
                    ));
                    planes.push(generate_plane(
                        PlaneType::VT,
                        &rect_lights,
                        id_emitter,
                        sampler,
                        m,
                    ));
                    planes.push(generate_plane(
                        PlaneType::UT,
                        &rect_lights,
                        id_emitter,
                        sampler,
                        m,
                    ));
                }
                SinglePlaneStrategy::UAlpha
                | SinglePlaneStrategy::ContinousMIS
                | SinglePlaneStrategy::SMISAll(_)
                | SinglePlaneStrategy::SMISJacobian(_)
                | SinglePlaneStrategy::ProxySample => {
                    planes.push(generate_plane(
                        PlaneType::UAlphaT,
                        &rect_lights,
                        id_emitter,
                        sampler,
                        m,
                    ));
                }
            }

            number_plane_gen += 1;
        }

        let s_size=20;
        let e_size=100;
        let x_size=1000;
        // build h_x cache with size as rect_lights.len()*s_size*8s_size*e_size*x_size
        let mut h_x = vec![vec![vec![vec![vec![0.0f32;x_size];e_size];8*s_size];s_size];rect_lights.len()];
        let mut b_max=vec![vec![vec![vec![0.0f32;e_size];8*s_size];s_size];rect_lights.len()];
        match self.strategy {
            SinglePlaneStrategy::ProxySample => {
                let light_length = rect_lights.len();
                
                for light_index in 0..light_length{
                    for u_index in 0..s_size as usize{
                        for v_index in 0..8*s_size as usize{
                            for e_index in 0..e_size{
                                let mut unnormalize_h_x=vec![0.0f32;x_size];
                                let mut accumulate=0.0f32;
                                for x_index in 0..x_size {
                                    //Now f(x)=Asin(pi*x)+B*cos(pi*x)
                                    let bias=PI/(e_size as f32)*(e_index as f32+0.5);
                                    let x=(x_index as f32+0.5)/(x_size as f32);
                                    let rect_light=&rect_lights[light_index];
                                    let dx: f32=(u_index as f32+0.5)/(s_size as f32);
                                    let dy: f32=(v_index as f32+0.5)/(8.0*s_size as f32);

                                    //calculate the length of the segment
                                    let plane2d_its = |d: Vector2<f32>, o: Point2<f32>| {
                                        let t_0 = (-o.to_vec()).div_element_wise(d);
                                        let t_1 = (Vector2::new(rect_light.u_l, rect_light.v_l) - o.to_vec()).div_element_wise(d);
                                        let t_max_coord = Vector2::new(t_0.x.max(t_1.x), t_0.y.max(t_1.y));
                                        o + d * t_max_coord.x.min(t_max_coord.y)
                                    };
                    
                                    let o_plane = Point2::new(dx*rect_light.u_l, dy*rect_light.v_l);
                                    let alpha = std::f32::consts::PI * x;
                                    let d_plane: Vector2<f32> = Vector2::new(alpha.cos(), alpha.sin());
                                    let p1_2d = plane2d_its(d_plane, o_plane);
                                    let p2_2d = plane2d_its(-d_plane, o_plane);
                                    let uu=f32::sin(PI*x+bias).abs()*(p1_2d-p2_2d).magnitude();
                                    unnormalize_h_x[x_index]=uu;
                                    accumulate+=uu;
                                }
                                b_max[light_index][u_index][v_index][e_index]=accumulate;
                                unnormalize_h_x.iter_mut().for_each(|x|{*x/=accumulate;});
                                h_x[light_index][u_index][v_index][e_index]=unnormalize_h_x;
                                //println!("h_x[{}][{}][{}][{}]={:?}",light_index,u_index,v_index,e_index,h_x[light_index][u_index][v_index][e_index]);
                            }
                        }
                    }
                }
            }
            _=> {}
        };
        // Build the BVH to speedup the computation...
        let bvh_plane = BHVAccel::create(planes);

        // Generate the image block to get VPL efficiently
        let buffernames = vec![String::from("primal")];
        let mut image_blocks = generate_img_blocks(scene, sampler, &buffernames);
        let sample_num=Mutex::new(0.0);
        // Gathering all planes
        info!("Gathering Single planes...");
        let progress_bar = Mutex::new(ProgressBar::new(image_blocks.len() as u64));
        let pool = generate_pool(scene);
        let phase_function = PhaseFunction::Isotropic();
        pool.install(|| {
            image_blocks
                .par_iter_mut()
                .for_each(|(im_block, sampler_ray)| {
                    // Sensor sampler is not used as to produce SMIS results, no AA was used.
                    let mut _sampler_ray = samplers::independent::IndependentSampler::from_seed(
                        (im_block.pos.x + im_block.pos.y) as u64,
                    );
                    let mut sampler_ecmis = samplers::independent::IndependentSampler::from_seed(
                        (im_block.pos.x + im_block.pos.y) as u64,
                    );
                    let mut block_sampel_num=0;
                    for ix in 0..im_block.size.x {
                        for iy in 0..im_block.size.y {
                            for _ in 0..scene.nb_samples {
                                let (ix_c, iy_c) = (ix + im_block.pos.x, iy + im_block.pos.y);
                                let pix = Point2::new(
                                    ix_c as f32 + sampler_ray.next(),
                                    iy_c as f32 + sampler_ray.next(),
                                );
                                let mut ray = scene.camera.generate(pix);

                                // Get the max distance
                                let max_dist = match accel.trace(&ray) {
                                    Some(x) => x.dist,
                                    None => std::f32::MAX,
                                };
                                ray.tfar = max_dist;

                                // Now gather all planes
                                let mut c = Color::value(0.0);
                                for (plane_its, b_id) in bvh_plane.gather(&ray) {
                                    let plane = &bvh_plane.elements[b_id];
                                    // This code is if we do not use BVH
                                    // for plane in &planes {
                                    // let plane_its = plane.intersection(&ray);
                                    // if plane_its.is_none() {
                                    // 	continue;
                                    // }
                                    // let plane_its = plane_its.unwrap();

                                let p_hit = ray.o + ray.d * plane_its.t_cam;
                                let p_light = plane
                                    .light_position(&rect_lights[plane.id_emitter], &plane_its);
                                let rect_light = &rect_lights[plane.id_emitter];
                                if accel.visible(&p_hit, &p_light) {
                                    let transmittance = {
                                        let mut ray_tr = Ray::new(ray.o, ray.d);
                                        ray_tr.tfar = plane_its.t_cam;
                                        m.transmittance(ray_tr)
                                    };
                                    let rho = phase_function
                                        .eval(&(-ray.d), &(p_light - p_hit).normalize());
                                    let w: f32 = match self.strategy {
                                        SinglePlaneStrategy::UT
                                        | SinglePlaneStrategy::UV
                                        | SinglePlaneStrategy::VT
                                        | SinglePlaneStrategy::UAlpha
                                        | SinglePlaneStrategy::ContinousMIS
                                        | SinglePlaneStrategy::SMISAll(_)
                                        | SinglePlaneStrategy::SMISJacobian(_)
                                        | SinglePlaneStrategy::ProxySample => 1.0,
                                        SinglePlaneStrategy::Average => 1.0 / 3.0,
                                        SinglePlaneStrategy::DiscreteMIS => {
                                            // Need to compute all possible shapes
                                            let d = p_hit - p_light;
                                            // TODO: Not used
                                            let t_sampled = d.magnitude();
                                            let d = d / t_sampled;
                                            let planes = [
                                                SinglePhotonPlane::new(
                                                    PlaneType::UV,
                                                    &rect_light,
                                                    d,
                                                    plane.sample,
                                                    0.0,
                                                    t_sampled,
                                                    plane.id_emitter,
                                                    m.sigma_s,
                                                ),
                                                SinglePhotonPlane::new(
                                                    PlaneType::UT,
                                                    &rect_light,
                                                    d,
                                                    plane.sample,
                                                    0.0,
                                                    t_sampled,
                                                    plane.id_emitter,
                                                    m.sigma_s,
                                                ),
                                                SinglePhotonPlane::new(
                                                    PlaneType::VT,
                                                    &rect_light,
                                                    d,
                                                    plane.sample,
                                                    0.0,
                                                    t_sampled,
                                                    plane.id_emitter,
                                                    m.sigma_s,
                                                ),
                                            ];
                                            // FIXME: Normally this code is unecessary
                                            // 	As we can reuse the plane retrived.
                                            // 	However, it seems to have a miss match between photon planes
                                            //	contribution calculation.
                                            let debug_id = match plane.plane_type {
                                                PlaneType::UV => 0,
                                                PlaneType::UT => 1,
                                                PlaneType::VT => 2,
                                                _ => unimplemented!(),
                                            };

                                            let w = planes[debug_id].contrib(&ray.d).avg().powi(-1)
                                                / planes
                                                    .iter()
                                                    .map(|p| {
                                                        let c = p.contrib(&ray.d).avg();
                                                        if c != 0.0 && c.is_finite() {
                                                            c.powi(-1)
                                                        } else {
                                                            0.0
                                                        }
                                                    })
                                                    .sum::<f32>();
                                            if w.is_finite() {
                                                w
                                            } else {
                                                0.0
                                            }
                                        }
                                    };

                                    // Compute the plane weighted constribution
                                    let contrib = match self.strategy {
                                        // Deng et al. CMIS
                                        SinglePlaneStrategy::ContinousMIS => {
                                            // Here we use their integration from
                                            // Normally, all the jacobian simplifies
                                            // So it is why we need to have a special estimator
                                            block_sampel_num+=1;
                                            1.0 / ((2.0 / std::f32::consts::PI)
                                                * (rect_light.u.cross(plane.d1).dot(ray.d).powi(2)
                                                    + rect_light.v.cross(plane.d1).dot(ray.d).powi(2))
                                                .sqrt())
                                                * plane.weight // The original jacobian cancel out
                                        }
                                        // SMIS
                                        SinglePlaneStrategy::SMISAll(n_samples)
                                        | SinglePlaneStrategy::SMISJacobian(n_samples) => {
                                            assert!(n_samples > 0);

                                            // Compute wrap random number for generating fake planes that generate same
                                            // path configuration. Indeed, we need to be sure that the new planes alpha plane
                                            // cross the same point on the light source
                                            let mut sample_wrap = {
                                                let p_l = p_light - rect_light.o;
                                                Point2::new(
                                                    p_l.dot(rect_light.u) / rect_light.u_l,
                                                    p_l.dot(rect_light.v) / rect_light.v_l,
                                                )
                                            };

                                            // Mitigate floating point precision issue
                                            sample_wrap.x = clamp(sample_wrap.x, 0.0, 1.0);
                                            sample_wrap.y = clamp(sample_wrap.y, 0.0, 1.0);

                                            // --------------------
                                            // Create N-1 fake planes with the same path configuration
                                            // using stratified sampling.

                                            // We will compute the inverse norm of the SMIS weight
                                            // We initialize the variable with the actual plane that
                                            // the ray intersected
                                            let mut inv_norm = match self.strategy {
                                                SinglePlaneStrategy::SMISAll(_v) => {
                                                    // J(..) * p(e)
                                                    // we do not include A (light source area)
                                                    // as it cancel out
                                                    plane.d1.cross(plane.d0).dot(ray.d).abs()
                                                        * plane.length0
                                                }
                                                SinglePlaneStrategy::SMISJacobian(_v) => {
                                                    // J(..)
                                                    plane.d1.cross(plane.d0).dot(ray.d).abs()
                                                }
                                                _ => panic!("Unimplemented"),
                                            };

                                            // Generate the other planes
                                            let offset = plane.sample_alpha;
                                            for i in 0..(n_samples - 1) {
                                                // The new alpha
                                                let new_alpha = {
                                                    if self.stratified {
                                                        (offset
                                                            + ((i + 1) as f32 / n_samples as f32))
                                                            % 1.0
                                                    } else {
                                                        sampler_ecmis.next()
                                                    }
                                                };
                                                assert!(new_alpha >= 0.0 && new_alpha <= 1.0);

                                                // This construct a new plane (call previous snipped)
                                                // to generate the plane
                                                let new_plane = SinglePhotonPlane::new(
                                                    PlaneType::UAlphaT,
                                                    &rect_light,
                                                    plane.d1,
                                                    sample_wrap, // Random number used to sample the point on the light source
                                                    new_alpha,
                                                    0.0,
                                                    plane.id_emitter,
                                                    m.sigma_s,
                                                );

                                                // Accumulate the weight
                                                inv_norm += match self.strategy {
                                                    SinglePlaneStrategy::SMISAll(_v) => {
                                                        new_plane
                                                            .d1
                                                            .cross(new_plane.d0)
                                                            .dot(ray.d)
                                                            .abs()
                                                            * new_plane.length0
                                                    }
                                                    SinglePlaneStrategy::SMISJacobian(_v) => {
                                                        new_plane
                                                            .d1
                                                            .cross(new_plane.d0)
                                                            .dot(ray.d)
                                                            .abs()
                                                    }
                                                    _ => panic!("Unimplemented"),
                                                };
                                            }

                                            // The SMIS weight and the plane contribution
                                            // First each weight have J(..) * p(e) or J(..) at the numerator
                                            // however this factor cancel out with the plane evaluation: f(..) / (J(..) * p(e))
                                            // With that, all planes contribute to the same :)
                                            let w = 1.0 / inv_norm;
                                            let contrib = match self.strategy {
                                                SinglePlaneStrategy::SMISAll(_v) => {
                                                    // Note we need to also cancel out lenght0
                                                    // as it will cancel out inside the SMIS_weight
                                                    plane.weight * plane.length0 * n_samples as f32
                                                }
                                                SinglePlaneStrategy::SMISJacobian(_v) => {
                                                    plane.weight * n_samples as f32
                                                }
                                                _ => panic!("Unimplemented"),
                                            };
                                            block_sampel_num+=n_samples;
                                            w * contrib
                                        }

                                        // Proxy sample
                                        SinglePlaneStrategy::ProxySample => {
                                            let mut sample_wrap = {
                                                let p_l = p_light - rect_light.o;
                                                Point2::new(
                                                    p_l.dot(rect_light.u) / rect_light.u_l,
                                                    p_l.dot(rect_light.v) / rect_light.v_l,
                                                )
                                            };

                                            // Mitigate floating point precision issue
                                            sample_wrap.x = clamp(sample_wrap.x, 0.0, 1.0);
                                            sample_wrap.y = clamp(sample_wrap.y, 0.0, 1.0);
                                            
                                            let mut inv_norm:f32=0.0;
                                            if self.stratified{
                                                let mut A=rect_light.v.cross(plane.d1).dot(ray.d);
                                                let mut B=rect_light.u.cross(plane.d1).dot(ray.d);
                                                if B<0.0 {
                                                    A=-A;
                                                    B=-B;
                                                }
                                                
                                                let dx: f32=rect_light.u_l*sample_wrap.x;
                                                let dy=rect_light.v_l*sample_wrap.y;
                                                let uu=direct_int(dx, rect_light.u_l-dx, rect_light.v_l-dy, A, B)
                                                +direct_int(rect_light.u_l-dx, dx, dy, A, B);
                                                inv_norm=1.0/uu;
                                                block_sampel_num+=1;
                                            }
                                            else 
                                            {
                                                let mut num_pos_stack = 1;
                                                let mut max_stack = 1000;

                                                //Now f(x)=Asin(pi*x)+B*cos(pi*x)
                                                let mut A=rect_light.v.cross(plane.d1).dot(ray.d);
                                                let mut B=rect_light.u.cross(plane.d1).dot(ray.d);
                                                if B<0.0 {
                                                    A=-A;
                                                    B=-B;
                                                }
                                                // let b=b_max[plane.id_emitter];
                                                let amplitude=(A.powf(2.0)+B.powf(2.0)).sqrt();
                                                let bias=f32::atan2(B,A);
                                                // let mut b=2.0*amplitude/PI*(rect_light.u_l.powf(2.0)+rect_light.v_l.powf(2.0)).sqrt();
                                                let a_index={
                                                    match (sample_wrap.x*(s_size as f32)) as usize{
                                                        u if u==s_size =>s_size-1,
                                                        u=>u,
                                                    }
                                                };
                                                let b_index={
                                                    match (sample_wrap.y*(8.0*s_size as f32)) as usize{
                                                        v if v==8*s_size =>8*s_size-1,
                                                        v=>v,
                                                    }
                                                };
                                                let e_index={
                                                    match (bias/PI*(e_size as f32)) as usize{
                                                        e if e==e_size =>e_size-1,
                                                        e=>e,
                                                    }
                                                };
                                                let alpha=0.0;//uniform_sample rate
                                                let b=amplitude*b_max[plane.id_emitter][a_index][b_index][e_index]/x_size as f32;
                                                //println!("b:{}",b);
                                                // let mut b=2.0*amplitude/PI*(rect_light.u_l.powf(2.0)+rect_light.v_l.powf(2.0)).sqrt();
                                                // let dx: f32=rect_light.u_l*sample_wrap.x;
                                                // let dy=rect_light.v_l*sample_wrap.y;
                                                // let uu=direct_int(dx, rect_light.u_l-dx, rect_light.v_l-dy, A, B)
                                                // +direct_int(rect_light.u_l-dx, dx, dy, A, B);
                                                //inv_norm=1.0/uu;

                                                
                                                //test code base
                                                // let mut v_c=vec![0.0f32;x_size];
                                                // let mut c_c=vec![0.0f32;x_size];
                                                // for n in 0..x_size{
                                                //     let x= (n as f32+0.5)/x_size as f32;
                                                //     let new_plane = SinglePhotonPlane::new(
                                                //         PlaneType::UAlphaT,
                                                //         &rect_light,
                                                //         plane.d1,
                                                //         sample_wrap,
                                                //         x,
                                                //         0.0,
                                                //         plane.id_emitter,
                                                //         m.sigma_s,
                                                //     );
                                                //     let ff=new_plane
                                                //     .d1
                                                //     .cross(new_plane.d0)
                                                //     .dot(ray.d).abs()*new_plane.length0;
                                                //     v_c[n]=ff;

                                                //     let c_sample_wrap=Point2::new(
                                                //         (a_index as f32+0.5)/s_size as f32,
                                                //         (b_index as f32+0.5)/(8.0*s_size as f32),
                                                //     );
                                                //     let new_plane = SinglePhotonPlane::new(
                                                //         PlaneType::UAlphaT,
                                                //         &rect_light,
                                                //         plane.d1,
                                                //         c_sample_wrap,
                                                //         x,
                                                //         0.0,
                                                //         plane.id_emitter,
                                                //         m.sigma_s,
                                                //     );
                                                //     let bias=PI/(e_size as f32)*(e_index as f32+0.5);
                                                //     let ff=new_plane.length0*amplitude*f32::sin(PI*x+bias).abs();
                                                //     c_c[n]=ff;
                                                // }
                                                // println!("v_c={:?}\nh_x={:?}\nc_c={:?}",v_c,h_x[plane.id_emitter][a_index][b_index][e_index],c_c);
                                                // exit(0);
                                                // let tao=0.5;
                                                // let normalization_factor={
                                                //     let normal = Normal::new(0.5, tao as f64).unwrap();
                                                //     (normal.cdf(1.0) - normal.cdf(0.0)) as f32
                                                // };
                                                // b=b*0.5+(1.0*PI).sqrt()*tao*rect_light.u_l.max(rect_light.v_l)/(-0.125/tao.powf(2.0)).exp();
                                                //println!("b:{} b1:{}",b,b1);
                                                //b=b.max(b1);
                                                
                                                let intergrate={
                                                    let u=rect_light.u_l;
                                                    let v=rect_light.v_l;
                                                    let a=rect_light.u_l*sample_wrap.x;
                                                    let b=rect_light.v_l*sample_wrap.y;
                                                    let rectangles = [(a, b), (u - a, b), (a, v - b), (u - a, v - b)];
                                                    let mut intergrate=0.0;
                                                    for &(m, n) in &rectangles {
                                                        intergrate+=calculate_f(m, n);
                                                        intergrate+=calculate_f(n, m);
                                                    }
                                                    intergrate/PI
                                                };  
                                                //println!("b:{}",b);
                                                b=2.0/(1.0/intergrate+1.0/b);
                                                // b=1.0;
                                                 while num_pos_stack != 0 && max_stack > 0 {
                                                    let mut sign = if num_pos_stack > 0 { 1.0 } else { -1.0 };
                                                    num_pos_stack -= sign as i32;

                                                    // let next_sample_normal=generate_normal_sample(tao);
                                                    // let plane_normal = SinglePhotonPlane::new(
                                                    //     PlaneType::UAlphaT,
                                                    //     &rect_light,
                                                    //     plane.d1,
                                                    //     sample_wrap, // Random number used to sample the point on the light source
                                                    //     next_sample_normal,
                                                    //     0.0,
                                                    //     plane.id_emitter,
                                                    //     m.sigma_s,
                                                    // );
                                                    // let f0: f32 = plane_normal
                                                    //             .d1
                                                    //             .cross(plane_normal.d0)
                                                    //             .dot(ray.d)
                                                    //             .abs()*plane_normal.length0;

                                                    // let next_sample_sin=generate_sin_sample(bias);
                                                    // let plane_sin=SinglePhotonPlane::new(
                                                    //     PlaneType::UAlphaT,
                                                    //     &rect_light,
                                                    //     plane.d1,
                                                    //     sample_wrap, // Random number used to sample the point on the light source
                                                    //     next_sample_sin,
                                                    //     0.0,
                                                    //     plane.id_emitter,
                                                    //     m.sigma_s,
                                                    // );
                                                    // let f1: f32 = plane_sin
                                                    //             .d1
                                                    //             .cross(plane_sin.d0)
                                                    //             .dot(ray.d)
                                                    //             .abs()*plane_sin.length0;

                                                    // let next_sample_len=generate_len_sample(rect_light.u_l, rect_light.v_l, rect_light.u_l*sample_wrap.x, rect_light.v_l*sample_wrap.y);
                                                    // let plane_len=SinglePhotonPlane::new(
                                                    //     PlaneType::UAlphaT,
                                                    //     &rect_light,
                                                    //     plane.d1,
                                                    //     sample_wrap, // Random number used to sample the point on the light source
                                                    //     next_sample_len,
                                                    //     0.0,
                                                    //     plane.id_emitter,
                                                    //     m.sigma_s,
                                                    // );
                                                    // let f2: f32 = plane_len
                                                    //             .d1
                                                    //             .cross(plane_len.d0)
                                                    //             .dot(ray.d)
                                                    //             .abs()*plane_len.length0;

                                                    let (next_sample_mis,p_mis)={
                                                        if sampler_ecmis.next()<alpha{
                                                            (sampler_ecmis.next(),1.0f32)
                                                        }
                                                        else{
                                                            let weights=(h_x[plane.id_emitter][a_index][b_index][e_index]).clone();
                                                            let dist = match WeightedIndex::new(&weights) {
                                                                Ok(d) => d,
                                                                Err(e) => {
                                                                    // 打印错误信息和权重数组
                                                                    eprintln!("Weights: {:?}", weights);
                                                                    exit(0);
                                                                }
                                                            };
                                                            
                                                            let mut rng = thread_rng();
                                                            let sampled_index = dist.sample(&mut rng);
                                                            ((sampled_index as f32+sampler_ecmis.next())/x_size as f32,h_x[plane.id_emitter][a_index][b_index][e_index][sampled_index]*x_size as f32)
                                                        }              
                                                    };
                                                    let plane_mis=SinglePhotonPlane::new(
                                                        PlaneType::UAlphaT,
                                                        &rect_light,
                                                        plane.d1,
                                                        sample_wrap, // Random number used to sample the point on the light source
                                                        next_sample_mis,
                                                        0.0,
                                                        plane.id_emitter,
                                                        m.sigma_s,
                                                    );
                                                    let f3: f32 = plane_mis
                                                                .d1
                                                                .cross(plane_mis.d0)
                                                                .dot(ray.d)
                                                                .abs()*plane_mis.length0;
                                                    let ff=f3/p_mis;
                                                    
                                                    // let p0_normal=normal_pdf(next_sample_normal,tao,normalization_factor);
                                                    // let p0_sin=sin_pdf(next_sample_normal,bias);
                                                    // let p1_normal=normal_pdf(next_sample_sin,tao,normalization_factor);
                                                    // let p1_sin=sin_pdf(next_sample_sin,bias);
                                                    // let p1_len=len_pdf(intergrate,plane_sin.length0);
                                                    // let p2_sin=sin_pdf(next_sample_len,bias);
                                                    // let p2_len=len_pdf(intergrate,plane_len.length0);

                                                    //let ff=f0/(p0_normal+p0_sin)+f1/(p1_normal+p1_sin);
                                                    // let ff=f1/(p1_sin+p1_len)+f2/(p2_sin+p2_len);
                                                    //let ff=f2/p2_len;

                                                    let g0=1.0-ff/b;
                                                    //assert!(g0>=0.0);
                                                    //println!("{} {}",g0,b);
                                                    inv_norm += 1.0 / b * sign; 
                                                    if g0<0.0 {
                                                        sign*=-1.0;
                                                    }
                                                    let r0 = g0.abs();
                                                    let mut r_round = r0.floor();
                                                    if sampler_ecmis.next() < r0 - r_round{
                                                        r_round += 1.0;
                                                    }
                                                    num_pos_stack += sign as i32 * r_round as i32;
                                                    max_stack -= 1;
                                                }
                                                // 
                                                //println!("Iter:{}",1000-max_stack);
                                                block_sampel_num+=1000-max_stack;
                                            }
                                            inv_norm*plane.weight*plane.length0
                                        }
                                        // Default: evaluate and weight the contrib
                                        // plane.contrib(..) =  plane.weight / jacobian
                                        // w: other MIS compute above
                                        _ => w * plane.contrib(&ray.d), // Do nothing and just evaluate the plane
                                    };

                                    // Compute the rest of the term
                                    // and accumulate them
                                    c += rho
                                        * transmittance
                                        * m.sigma_s
                                        * contrib
                                        * (emitters.len() as f32)
                                        * (1.0 / number_plane_gen as f32);
                                }
                            }
                                im_block.accumulate(
                                    Point2 { x: ix, y: iy },
                                    c,
                                    &"primal".to_owned(),
                                );
                            }
                        }
                    } // Image block
                    {
                        let mut num =sample_num.lock().unwrap();
                        *num+=block_sampel_num as f32;
                    }
                    im_block.scale(1.0 / (scene.nb_samples as f32));
                    {
                        progress_bar.lock().unwrap().inc();
                    }
                });
        });
        println!("\n Sample number_all:{}",sample_num.lock().unwrap());
        println!("Sample number:{}",*sample_num.lock().unwrap()/(scene.camera.size().x*scene.camera.size().y) as f32);
        let mut image =
            BufferCollection::new(Point2::new(0, 0), *scene.camera.size(), &buffernames);
        for (im_block, _) in &image_blocks {
            image.accumulate_bitmap(im_block);
        }
        image
    }
}
