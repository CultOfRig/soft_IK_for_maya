/*
Copyright 2019 Raffaele Fragapane

http://cultofrig.com
https://github.com/CultOfRig/soft_IK_for_maya

Redistribution and use in source and binary forms, with or without modification,
  are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
     this list of conditions and the following disclaimer
     in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
     may be used to endorse or promote products
     derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
   THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
   IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
   HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
   EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#pragma once

#include <vector>
#include <algorithm>

#ifdef _WIN64
#define _USE_MATH_DEFINES
#endif
#include <cmath>

#include <maya/MPxNode.h>

#include <maya/MFnNumericAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnEnumAttribute.h>

#include <maya/MFnNumericData.h>

#include <maya/MFloatVector.h>
#include <maya/MFloatMatrix.h>

static const int k_vchain_soft_id = 0x0012a040;
static const char* k_vchain_soft_name = "cor_vChain_multi_solve";


static const uint32_t k_output_count = 7;
const MFloatMatrix k_mtx44f_id = MFloatMatrix();


static const short k_solver_rigid = 0;
static const short k_solver_height_reduction = 1;
static const short k_solver_height_lock = 2;
static const short k_solver_length_dampening= 3;


struct Triangle
{
    float a = 0.f;
    float b = 0.f;
    float c = 0.f;

    double alpha = 0.0;
    double beta = 0.0;
    double gamma = 0.0;

    float sin_alpha = 1.f;
    float cos_alpha = 0.f;

    float sin_beta = 0.f;
    float cos_beta = 1.f;

    float sin_gamma = 0.5f;
    float cos_gamma = 0.5f;
};


class VChain_soft : public MPxNode
{
public:
    VChain_soft();

    virtual SchedulingType schedulingType()const
    {
        return SchedulingType::kParallel;
    }

    static void* creator(){return new VChain_soft();};

    static MStatus initialize();

    virtual MStatus compute(const MPlug& plug, MDataBlock& data_block);

    //  INPUT
    static MObject in_parent_frame;

    static MObject in_root;
    static MObject in_handle;

    static MObject in_up_vector;
    static MObject in_roll_angle;

    static MObject in_length1;
    static MObject in_length2;

    static MObject in_max_compression;
    static MObject in_min_extension;
    static MObject in_max_extension;

    static MObject in_solver_choice;

    static MObject in_hierarchical_mode;



    // OUTPUT
    static MObject out_result;

private:
    static const uint32_t output_count = k_output_count;
    std::vector<MFloatMatrix> mtx44fa_result;

    Triangle rigid_tri;
    Triangle soft_tri;
};

