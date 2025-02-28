//
// Created by dhruv on 5/26/23.
//

#include <iostream>
#include <DendriteUtils.h>
#include <point.h>
#include <TalyEquation.h>
#include <DeepTraceInputData.h>
#include <IO/VTU.h>
#include "SDARefine.h"
#include "talyfem/utils/timers.h"
#include <onnxruntime_cxx_api.h>
#include <torch/script.h>

void performNoChangeRefinement(DA *& octDA, DistTREE & distTree, DomainExtents &domainInfo,
                               DENDRITE_UINT maxLevel, SubDomain & subDomain){
    // Use this function to perform no refinement but to remove the physical boundary octants in case of
    // retain inside case
    SDARefine refine(octDA, distTree.getTreePartFiltered(),domainInfo,maxLevel);
    DA * newDA = refine.getForceRefineSubDA(distTree,0.01);
    std::swap(octDA, newDA);
    delete newDA;
    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainInfo);

}

void performRefinement(DA *& octDA, DistTREE & distTree, DomainExtents &domainInfo,
                               DENDRITE_UINT maxLevel, SubDomain & subDomain){
    while(true) {
        SDARefine refine(octDA, distTree.getTreePartFiltered(), domainInfo, maxLevel);

        DA *newDA = refine.getRefineSubDA(distTree, 0.01,RefinementStrategy::FULL_DOMAIN,ot::RemeshPartition::SurrogateInByOut);
        if(newDA == nullptr){
            break;
        }
        std::swap(octDA, newDA);
        delete newDA;
        subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainInfo);
    }

}

///

// Function to generate random points in the range [minValue, maxValue]
float generateRandomFloat(float minValue, float maxValue) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(minValue, maxValue);
    return dis(gen);
}

// Function to generate random 3D points
std::vector<float> generateRandomPoints(int numPoints) {
    std::vector<float> points;

    // Generate random x, y, z coordinates for each point
    for (int i = 0; i < numPoints; ++i) {
        points.push_back(generateRandomFloat(0.0f, 1.0f));  // Random x coordinate
        points.push_back(generateRandomFloat(0.0f, 1.0f));  // Random y coordinate
        points.push_back(generateRandomFloat(0.0f, 1.0f));  // Random z coordinate
    }

    return points;
}


///
int main(int argc, char *argv[]) {
  dendrite_init(argc, argv);
  /// read parameters from config.txt
  DeepTraceInputData idata;

  /// use config.txt first
  std::ifstream configFile("config.txt");
  DENDRITE_UINT eleOrder = 0;
  DENDRITE_UINT uniLevel = 0;
  DENDRITE_UINT bdLevel = 0;
  if (configFile.good()) {
    if (!idata.ReadFromFile()) {  /// read from file named "config.txt"
      throw std::runtime_error("[ERR] Error reading input data, check the config file!");
    }
    if (!idata.CheckInputData()) {
      throw std::runtime_error("[ERR] Problem with input data, check the config file!");
    }
    eleOrder = idata.basisFunction;
    uniLevel = idata.refine_lvl_uni;
    bdLevel = idata.refine_lvl_bd;
  } else {
    TALYFEMLIB::PrintStatus("Fix the config file!");
    return -1;
  }

  TALYFEMLIB::PrintInfo("Total number of processor = ", TALYFEMLIB::GetMPISize());
  TALYFEMLIB::PrintInfo("size of DendroInt ", sizeof(DendroIntL));
  TALYFEMLIB::PrintInfo("size of PetscInt ", sizeof(PetscInt));

  TALYFEMLIB::PrintStatus("eleOrder ", eleOrder);
  TALYFEMLIB::PrintStatus("Level ", uniLevel);
  TALYFEMLIB::PrintStatus("DIM =  ", DIM);


    /// Physical dimensions.
    DomainInfo cubeDomain, physDomain;
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
        cubeDomain.min[dim] = idata.cubeDomainMin[dim];
        cubeDomain.max[dim] = idata.cubeDomainMax[dim];
        physDomain.min[dim] = idata.physDomainMin[dim];
        physDomain.max[dim] = idata.physDomainMax[dim];
    }
    DomainExtents domainExtents(cubeDomain,physDomain);
    SubDomain subDomain(domainExtents);

        MPITimer timer("Run Time");
        timer.Start();







#if(DIM==2)
    {
        VOXEL::Circle circle(idata.circleCenter, idata.circleRadius, idata.circleRetainSide);
        CreateMesh(idata, circle, uniLevel, bdLevel, eleOrder, "circle");
    }

    {
        VOXEL::Box box(idata.rectLowerLeft, idata.rectUpperRight, idata.rectRetainSide);
        CreateMesh(idata,box,uniLevel, bdLevel, eleOrder,"rectangle");
    }
    {
        GEOMETRY::MSH msh(idata.mshFileName.c_str());
        Point<DIM> temp(0.0);
        GEOMETRY::Geometry *mshGeo = new GEOMETRY::Geometry(&msh, temp, idata.mshRetainSide);
        CreateMesh(idata,mshGeo,uniLevel, bdLevel, eleOrder,"msh");
        delete mshGeo;
    }

#endif
#if(DIM==3)
//    {
//        VOXEL::Sphere sphere(idata.sphereCenter, idata.sphereRadius, idata.sphereRetainSide);
//        CreateMesh(idata,sphere,uniLevel, bdLevel, eleOrder,"sphere");
//    }
//
//    {
//        VOXEL::Box box(idata.boxLowerLeft, idata.boxUpperRight, idata.boxRetainSide);
//        CreateMesh(idata,box,uniLevel, bdLevel, eleOrder,"box");
//    }







//    // Call the function you want to measure
//
//    // Initialize the ONNX Runtime environment
//    Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "test");
//
//    // Set up session options
//    Ort::SessionOptions session_options;
//    session_options.SetIntraOpNumThreads(1);
//
//
//
//    // Create a session with the loaded model
//    Ort::Session session(env, idata.ModelFileName, session_options);
//
//    // Define names of the input and output tensors
//    const char* input_name = "onnx::Gemm_0"; // Replace with the actual input name obtained from the model
//    const char* output_name = "46"; // Replace with the actual output name obtained from the model





//        // Deserialize the ScriptModule from a file using torch::jit::load().
//        torch::jit::script::Module module = torch::jit::load("model.pt");
//
//        // Create a tensor with the specified coordinates (0.5, 0.5, 0.5).
//        torch::Tensor coordinates = torch::full({1, 3}, 0.5);
//
//        // Wrap the input tensor in a vector of IValue.
//        std::vector<torch::jit::IValue> inputs;
//        inputs.push_back(coordinates);
//
//        // Execute the model and turn its output into a tensor.
//        at::Tensor output = module.forward(inputs).toTensor();
//
//        // Print the output tensor.
//        std::cout << "Output Tensor:\n" << output << '\n';
//
//        // Just testing with the model in the directory for verification
//        at::Tensor output_python = at::tensor(0.1233);
//        std::cout << "Output Tensor from Python:\n" << output_python << '\n';


        const auto geomDef = [&](const double * coords) {
        // Convert the physical coordinates to image coordinates using idata.physDomainMin, idata.physDomainMax
        // and idata.image_size

//        MPITimer timer("DeepTrace");
//        timer.Start();
//
//        float x = coords[0];
//        float y = coords[1];
//        float z = coords[2];
//
//
//
//        // Prepare input tensor
//        std::vector<float> input_tensor_values = {x, y, z};
//        std::vector<int64_t> input_tensor_shape = {1, 3}; // Assuming the model expects a [1,3] shape tensor
//
//        // Create memory info to describe the allocation info
//        Ort::MemoryInfo memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
//
//        // Create input tensor
//        Ort::Value input_tensor = Ort::Value::CreateTensor<float>(
//                memory_info, input_tensor_values.data(), input_tensor_values.size(), input_tensor_shape.data(), input_tensor_shape.size());
//
////        // Check if the model expects the input tensor
////        if (!input_tensor.IsTensor()) {
////            std::cerr << "Input is not a tensor\n";
////            return -1;
////        }
//
//
//        // Perform the inference
//        auto output_tensors = session.Run(Ort::RunOptions{nullptr}, &input_name, &input_tensor, 1, &output_name, 1);
//
//
//        float* floatarr = output_tensors.front().GetTensorMutableData<float>();
////        std::cout << "Inference result: " << floatarr[0] << std::endl;
//
//        float dist = 0;
//
//        timer.Stop();
//        timer.PrintTotalTimeSeconds();
//
//        if (idata.stlRetainSide == RetainSide::IN){
//            dist = floatarr[0];
//        }
//        else{
//            dist = -floatarr[0];
//        }
//
//
//        if (dist > 0){
//
//            return ibm::Partition::IN;
//        }
        return ibm::Partition::OUT;


    };

        GEOMETRY::STL stl(idata.stlFileName.c_str());
        Point<DIM> temp2(0.0);
        GEOMETRY::Geometry *stlGeo = new GEOMETRY::Geometry(&stl, temp2, idata.stlRetainSide);



#endif

    if (idata.useDeepLearning) {
        subDomain.addObject(geomDef);
    }
    else{
        subDomain.addObject(stlGeo);
    }


    // Retain Function
    auto functionToRetain = [&](const double *physCoords,double physSize) {
        return (subDomain.functionToRetain(physCoords, physSize));
    };

    /// Create DA
    DistTREE distTree;
    DA * octDA = createSubDA(distTree, functionToRetain, uniLevel, eleOrder, 0.3);
    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainExtents);
    performRefinement(octDA, distTree, domainExtents, bdLevel, subDomain);
    performNoChangeRefinement(octDA, distTree, domainExtents, 5, subDomain);


    timer.Stop();
    timer.PrintTotalTimeSeconds();

    IO::writeBoundaryElements(octDA, distTree.getTreePartFiltered(), "output", "output", domainExtents);


    dendrite_finalize();
}