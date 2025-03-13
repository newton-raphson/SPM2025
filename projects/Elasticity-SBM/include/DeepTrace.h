//
// Created by samundra on 3/21/24.
//

#ifndef NS_LD_CONS_IMGA_DEEPTRACE_DEEPTRACE_H
#define NS_LD_CONS_IMGA_DEEPTRACE_DEEPTRACE_H

#include <onnxruntime_cxx_api.h>
#include <cmath>
#include <iostream>
#include <mutex>
#include <memory>


class DeepTrace {
public:
    DeepTrace(const std::string& input_name, const std::string& output_name, const std::string& model_path, double scale);
    DeepTrace(const std::string& model_path, double scale);

    // Delete copy constructor (Ort::Session is not copyable)
    DeepTrace(const DeepTrace& other) = delete;
    DeepTrace& operator=(const DeepTrace& other) = delete;

    double compute_signed_distance(const double* coords);
    void compute_distance_vector(const double* coords, double* distance_vector);

private:
    static std::shared_ptr<Ort::Session> session;  // Shared ONNX session
    static std::mutex session_mutex;  // Thread safety

    Ort::Env env;
    Ort::SessionOptions session_options;
    Ort::MemoryInfo memory_info;
    double scale;
    std::string input_name;
    std::string output_name;
    std::string model_path;
};

// Initialize static members
std::shared_ptr<Ort::Session> DeepTrace::session = nullptr;
std::mutex DeepTrace::session_mutex;

DeepTrace::DeepTrace(const std::string& input_name, const std::string& output_name, const std::string& model_path, double scale)
        : memory_info(Ort::MemoryInfo::CreateCpu(OrtDeviceAllocator, OrtMemTypeCPU)),
          env(ORT_LOGGING_LEVEL_WARNING, "DeepTrace"),
          scale(scale),
          input_name(input_name),
          output_name(output_name),
          model_path(model_path) {

    std::lock_guard<std::mutex> lock(session_mutex);  // Ensure thread safety

    if (!session) {  // Load model only once
        session_options.SetIntraOpNumThreads(1);
        session_options.SetInterOpNumThreads(1);
        session_options.DisablePerSessionThreads();
        session = std::make_shared<Ort::Session>(env, model_path.c_str(), session_options);
    }
}

DeepTrace::DeepTrace(const std::string& model_path, double scale)
        : DeepTrace("input", "output", model_path, scale) {}  // Delegate to main constructor

double DeepTrace::compute_signed_distance(const double *coords) {


//    scale the coordinates

//    x = (x - 0.5 )*0.85/0.167;
//    y = (y - 0.5 )*0.85/0.167;
//    std::vector<float> float_coords{static_cast<float>((coords[0]-0.5)*0.85*6), static_cast<float>((coords[1]-0.5)*0.85*6), static_cast<float>(coords[2])};
    std::vector<float> float_coords{static_cast<float>((coords[0]*scale)), static_cast<float>((coords[1]*scale)), static_cast<float>(coords[2]*scale)};
    // Prepare input tensor
//    std::vector<float> input_tensor_values{float_coords[0], float_coords[1], float_coords[2]};
    std::vector<int64_t> input_tensor_shape{1, 3}; // Assuming the model expects a [1,3] shape tensor
    Ort::MemoryInfo memory_info_c(Ort::MemoryInfo::CreateCpu(OrtDeviceAllocator, OrtMemTypeCPU));
    Ort::Value input_tensor = Ort::Value::CreateTensor<float>(
                memory_info_c, float_coords.data(), float_coords.size(), input_tensor_shape.data(), input_tensor_shape.size());
//        std::cout<<"\n";
//        std::cout<<input_tensor;
//        std::cout<<"\n";
    // Perform the inference
//    auto input_name = "input";
//    auto output_name="output";

    // Run inference
    const char* input_names[] = {input_name.c_str()};
    const char* output_names[] = {output_name.c_str()};
    auto output_tensors = session->Run(
            Ort::RunOptions{nullptr},
            input_names, &input_tensor, 1,
            output_names, 1
    );


    // Check if the inference succeeded
    if (output_tensors.empty()) {
        std::cerr << "Inference failed\n";
        return -1.0; // Or any other appropriate error handling
    }

    // Get the signed distance value from the output tensor
    float dist_float = *output_tensors.front().GetTensorData<float>();

    auto dist = static_cast<double>(dist_float); // Recast to double
    return dist/scale;

}

void DeepTrace::compute_distance_vector(const double *coords, double *gradient) {
    constexpr double epsilon = 1e-6; // A slightly larger epsilon for more accuracy with double precision
    constexpr double epsilon_2 = epsilon * epsilon;

    // Compute the signed distance for the given coordinates
    double base_distance = DeepTrace::compute_signed_distance(coords);

    double magnitude_squared = 0.0; // To compute the magnitude of the gradient

    for (int i = 0; i < DIM; ++i) {
        double perturbed_coords[DIM];
        std::copy(coords, coords + DIM, perturbed_coords); // Make a copy of the original coordinates

        // Perturb the i-th coordinate positively and compute signed distance
        perturbed_coords[i] += epsilon;
        double f_x_plus_h = compute_signed_distance(perturbed_coords);

        // Perturb the i-th coordinate negatively and compute signed distance
        perturbed_coords[i] -= 2 * epsilon; // Now, perturb in the negative direction
        double f_x_minus_h = compute_signed_distance(perturbed_coords);

        // Compute the numerical gradient using a more accurate central difference method
        double gradient_value = (8.0 * (f_x_plus_h - f_x_minus_h) - (compute_signed_distance(coords) - compute_signed_distance(coords))) / (12.0 * epsilon);

        gradient[i] = gradient_value;
        magnitude_squared += gradient_value * gradient_value;
    }
// if base_distance is positive flip the sign of the gradient
    if (base_distance > 0) {
        for (int i = 0; i < DIM; ++i) {
            gradient[i] = -gradient[i];
        }
    }
    // Compute the magnitude of the gradient
    double magnitude = std::sqrt(magnitude_squared);

    // Avoid division by zero, if magnitude is non-zero normalize it
    if (magnitude > 1e-10) {
        double scale_factor = base_distance / magnitude;
        for (int i = 0; i < DIM; ++i) {
            gradient[i] *= scale_factor;
        }
    }
}




//
//class DeepTrace
//{
//public:
//    DeepTrace(const char* input_name,const char* output_name,const char* model_path,double scale);
//    DeepTrace(const char* model_path, double scale);
//    DeepTrace(const DeepTrace &other); // Copy constructor
//    double compute_signed_distance(const double * coords);
//    void compute_distance_vector(const double * coords,double *distance_vector);
//
//
//private:
//    Ort::Env env;
//    Ort::SessionOptions session_options;
//    Ort::MemoryInfo memory_info;
//    Ort::Session session;
//    const char* input_name;
//    const char* output_name;
//    const char* model_path;
//    double scale;
//
//};
//
//DeepTrace::DeepTrace(const char *input_name, const char *output_name, const char *model_path,double scale)
//        : memory_info(Ort::MemoryInfo::CreateCpu(OrtDeviceAllocator, OrtMemTypeCPU)),
//          env(ORT_LOGGING_LEVEL_WARNING, "test"),
//          session(env,model_path , session_options),scale(scale)
//{
//
//
//    DeepTrace::input_name = input_name;
//    DeepTrace::output_name = output_name;
//    DeepTrace::model_path = model_path;
//    session_options.SetIntraOpNumThreads(1);
//
//
//}
//DeepTrace::DeepTrace( const char *model_path,double scale)
//        : memory_info(Ort::MemoryInfo::CreateCpu(OrtDeviceAllocator, OrtMemTypeCPU)),
//          env(ORT_LOGGING_LEVEL_WARNING, "test"),
//          session(env,model_path , session_options),scale(scale)
//{
//    DeepTrace::input_name = "input";
//    DeepTrace::output_name = "output";
//    DeepTrace::model_path = model_path;
//    session_options.SetIntraOpNumThreads(2);
//
//}
//double DeepTrace::compute_signed_distance(const double *coords) {
//
//
////    scale the coordinates
//
////    x = (x - 0.5 )*0.85/0.167;
////    y = (y - 0.5 )*0.85/0.167;
////    std::vector<float> float_coords{static_cast<float>((coords[0]-0.5)*0.85*6), static_cast<float>((coords[1]-0.5)*0.85*6), static_cast<float>(coords[2])};
//    std::vector<float> float_coords{static_cast<float>((coords[0]*scale)), static_cast<float>((coords[1]*scale)), static_cast<float>(coords[2]*scale)};
//    // Prepare input tensor
////    std::vector<float> input_tensor_values{float_coords[0], float_coords[1], float_coords[2]};
//    std::vector<int64_t> input_tensor_shape{1, 3}; // Assuming the model expects a [1,3] shape tensor
//    Ort::MemoryInfo memory_info_c(Ort::MemoryInfo::CreateCpu(OrtDeviceAllocator, OrtMemTypeCPU));
//        Ort::Value input_tensor = Ort::Value::CreateTensor<float>(
//                memory_info_c, float_coords.data(), float_coords.size(), input_tensor_shape.data(), input_tensor_shape.size());
////        std::cout<<"\n";
////        std::cout<<input_tensor;
////        std::cout<<"\n";
//    // Perform the inference
////    auto input_name = "input";
////    auto output_name="output";
//    auto output_tensors = session.Run(Ort::RunOptions{nullptr}, &input_name, &input_tensor, 1, &output_name, 1);
//
//    // Check if the inference succeeded
//    if (output_tensors.empty()) {
//        std::cerr << "Inference failed\n";
//        return -1.0; // Or any other appropriate error handling
//    }
//
//    // Get the signed distance value from the output tensor
//    float dist_float = *output_tensors.front().GetTensorData<float>();
//
//    auto dist = static_cast<double>(dist_float); // Recast to double
////    std::cout<<"Signed_Distance\n";
////    std::cout<<dist<<std::endl;
//    return dist/scale;
//
//}
//
//void DeepTrace::compute_distance_vector(const double *coords, double *gradient) {
//    constexpr double epsilon = 1e-6; // A slightly larger epsilon for more accuracy with double precision
//    constexpr double epsilon_2 = epsilon * epsilon;
//
//    // Compute the signed distance for the given coordinates
//    double base_distance = DeepTrace::compute_signed_distance(coords);
//
//    double magnitude_squared = 0.0; // To compute the magnitude of the gradient
//
//    for (int i = 0; i < DIM; ++i) {
//        double perturbed_coords[DIM];
//        std::copy(coords, coords + DIM, perturbed_coords); // Make a copy of the original coordinates
//
//        // Perturb the i-th coordinate positively and compute signed distance
//        perturbed_coords[i] += epsilon;
//        double f_x_plus_h = compute_signed_distance(perturbed_coords);
//
//        // Perturb the i-th coordinate negatively and compute signed distance
//        perturbed_coords[i] -= 2 * epsilon; // Now, perturb in the negative direction
//        double f_x_minus_h = compute_signed_distance(perturbed_coords);
//
//        // Compute the numerical gradient using a more accurate central difference method
//        double gradient_value = (8.0 * (f_x_plus_h - f_x_minus_h) - (compute_signed_distance(coords) - compute_signed_distance(coords))) / (12.0 * epsilon);
//
//        gradient[i] = gradient_value;
//        magnitude_squared += gradient_value * gradient_value;
//    }
//// if base_distance is positive flip the sign of the gradient
//    if (base_distance > 0) {
//        for (int i = 0; i < DIM; ++i) {
//            gradient[i] = -gradient[i];
//        }
//    }
//    // Compute the magnitude of the gradient
//    double magnitude = std::sqrt(magnitude_squared);
//
//    // Avoid division by zero, if magnitude is non-zero normalize it
//    if (magnitude > 1e-10) {
//        double scale_factor = base_distance / magnitude;
//        for (int i = 0; i < DIM; ++i) {
//            gradient[i] *= scale_factor;
//        }
//    }
//}
//
//
//



#endif //NS_LD_CONS_IMGA_DEEPTRACE_DEEPTRACE_H
