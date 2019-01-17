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


#include "include/vChain_nodes.h"

VChain_soft::VChain_soft()
{
    mtx44fa_result.resize(k_output_count, k_mtx44f_id);
}

MObject VChain_soft::in_parent_frame;

MObject VChain_soft::in_root;
MObject VChain_soft::in_handle;
MObject VChain_soft::in_up_vector;

MObject VChain_soft::in_roll_angle;

MObject VChain_soft::in_length1;
MObject VChain_soft::in_length2;

MObject VChain_soft::in_max_compression;
MObject VChain_soft::in_min_extension;
MObject VChain_soft::in_max_extension;

MObject VChain_soft::in_solver_choice;

MObject VChain_soft::in_hierarchical_mode;

MObject VChain_soft::out_result;


MStatus VChain_soft::initialize()
{
    MStatus status = MStatus::kFailure;

    MFnNumericAttribute n_attr;
    MFnUnitAttribute u_attr;
    MFnMatrixAttribute m_attr;
    MFnEnumAttribute e_attr;

    // output
    out_result = m_attr.create("result", "result",
                               MFnMatrixAttribute::kFloat,
                               &status);
    m_attr.setWritable(false);
    m_attr.setArray(true);
    m_attr.usesArrayDataBuilder(false);

    addAttribute(out_result);


    // input
    //    handles
    in_parent_frame = m_attr.create("parent_frame", "parent_frame",
                                    MFnMatrixAttribute::kFloat,
                                    &status);
    m_attr.setWritable(true);
    
    in_root = m_attr.create("root", "root",
                            MFnMatrixAttribute::kFloat,
                            &status);
    m_attr.setWritable(true);

    in_handle = m_attr.create("handle", "handle",
                               MFnMatrixAttribute::kFloat,
                               &status);
    m_attr.setWritable(true);

    in_up_vector = m_attr.create("up_vector", "up_vector",
                                 MFnMatrixAttribute::kFloat,
                                 &status);
    m_attr.setWritable(true);

    // lengths
    in_length1 = n_attr.create("length1", "length1",
                               MFnNumericData::kFloat, 2.f,
                               &status);
    n_attr.setSoftMin(.5f); n_attr.setSoftMax(10.f);
    n_attr.setMin(.05f);
    n_attr.setKeyable(true);
    n_attr.setWritable(true);

    in_length2 = n_attr.create("length2", "length2",
                               MFnNumericData::kFloat, 2.f,
                               &status);
    n_attr.setSoftMin(.5f); n_attr.setSoftMax(10.f);
    n_attr.setMin(.05f);
    n_attr.setKeyable(true);
    n_attr.setWritable(true);

    in_max_compression = n_attr.create("max_compression_ratio", "max_compression_ratio",
                                       MFnNumericData::kFloat, 0.f,
                                       &status);
    n_attr.setSoftMin(.05f);
    n_attr.setMin(.0f); n_attr.setMax(1.f);
    n_attr.setKeyable(true);
    n_attr.setWritable(true);

    in_min_extension = n_attr.create("min_extension_ratio", "min_extension_ratio",
                                     MFnNumericData::kFloat, 0.9f,
                                     &status);
    n_attr.setSoftMin(.8f);
    n_attr.setMin(.0f); n_attr.setMax(1.f);
    n_attr.setKeyable(true);
    n_attr.setWritable(true);

    in_max_extension = n_attr.create("max_extension_ratio", "max_extension_ratio",
                                     MFnNumericData::kFloat, 1.1f,
                                     &status);
    n_attr.setSoftMax(2.f);
    n_attr.setMin(1.f);
    n_attr.setKeyable(true);
    n_attr.setWritable(true);

    in_solver_choice = e_attr.create("solver_choice", "solver_choice", 0, &status);
    e_attr.setKeyable(true);
    e_attr.addField("rigid", k_solver_rigid);
    e_attr.addField("height_reduction", k_solver_height_reduction);
    e_attr.addField("height_lock", k_solver_height_lock);
    e_attr.addField("length_dampening", k_solver_length_dampening);

    in_hierarchical_mode = n_attr.create("hierarchical_output", "hierarchical_output",
                                         MFnNumericData::kBoolean, true, &status);
    n_attr.setKeyable(true);
    n_attr.setWritable(true);


    addAttribute(in_parent_frame);
    attributeAffects(in_parent_frame, out_result);

    addAttribute(in_root);
    attributeAffects(in_root, out_result);

    addAttribute(in_handle);
    attributeAffects(in_handle, out_result);

    addAttribute(in_up_vector);
    attributeAffects(in_up_vector, out_result);


    addAttribute(in_length1);
    attributeAffects(in_length1, out_result);

    addAttribute(in_length2);
    attributeAffects(in_length2, out_result);


    addAttribute(in_max_compression);
    attributeAffects(in_max_compression, out_result);

    addAttribute(in_min_extension);
    attributeAffects(in_min_extension, out_result);

    addAttribute(in_max_extension);
    attributeAffects(in_max_extension, out_result);


    addAttribute(in_solver_choice);
    attributeAffects(in_solver_choice, out_result);

    addAttribute(in_hierarchical_mode);
    attributeAffects(in_hierarchical_mode, out_result);

    return MStatus::kSuccess;
}


/// \brief Local use composition of row vectors into an affine matrix compatible with
///          Maya's expectations
///
/// \param x   tri-vector of floats for X.x, X.y, X.z row. 4th element of the row implied as 0 <br>
///            Intentionally chosen to be passed by value and not ref/ptr as it barely outsizes
///              the reference size
/// \param y   tri-vector of floats for Y.x, Y.y, Y.z row. 4th element of the row implied as 0 <br>
///            Intentionally chosen to be passed by value and not ref/ptr as it barely outsizes
///              the reference size
/// \param z   tri-vector of floats for Z.x, Z.y, Z.z row. 4th element of the row implied as 0 <br>
///            Intentionally chosen to be passed by value and not ref/ptr as it barely outsizes
///              the reference size
/// \param pos tri-vector of floats for D.x, D.y, D.z row. 4th element of the row implied as 1 <br>
///            Intentionally chosen to be passed by value and not ref/ptr as it barely outsizes
///              the reference size
/// \return    Maya compatible affine matrix of the transformation
///              aligned to the axes and position passed in the arguments
static
MFloatMatrix matrix_from_vectors(MFloatVector x,
                                 MFloatVector y, 
                                 MFloatVector z,
                                 MFloatVector pos = MFloatVector::zero)
{
    MFloatMatrix ret_mtx(k_mtx44f_id);

    // The following might be non ideal compared to memcpy
    //   but should have better pedagogic value and clarity for the reader
    ret_mtx[0][0] = x.x; ret_mtx[0][1] = x.y; ret_mtx[0][2] = x.z;
    ret_mtx[1][0] = y.x; ret_mtx[1][1] = y.y; ret_mtx[1][2] = y.z;
    ret_mtx[2][0] = z.x; ret_mtx[2][1] = z.y; ret_mtx[2][2] = z.z;
    ret_mtx[3][0] = pos.x; ret_mtx[3][1] = pos.y; ret_mtx[3][2] = pos.z; ret_mtx[3][3] = 1.f;

    return ret_mtx;
}


/// \brief
///
/// \param solution a pointer to the triangle data containing the solution to convert to matrices
/// \param mtx_basis Matrix containing the basis (atV and upV transform) to transform to
/// \param mtx_handle Matrix containing world space of the handle to re-align the last transform to
/// \param mtx44fa_result A sequence of matrices representing transforms
///                         at key points in the triangle and algined to the sides.<br>
///                       This is an output argument.
/// \return This is a state-modifying (side-effect only) function
static
void result_from_triangle_hierarchical(const Triangle* const solution,
                                       const MFloatMatrix& mtx_basis,
                                       const MFloatMatrix& mtx_handle,
                                       std::vector<MFloatMatrix>& mtx44fa_result) // output argument
{
    { // zeroth
      // this is going to be the plain IK basis
        mtx44fa_result[0] = mtx_basis;
    }

    { // first joint
        MFloatMatrix mtx_local(k_mtx44f_id);
        mtx_local[0][0] = solution->cos_beta; mtx_local[0][1] = solution->sin_beta;
        mtx_local[1][0] = -solution->sin_beta; mtx_local[1][1] = solution->cos_beta;

        mtx44fa_result[1] = mtx_local;
    }

    { // second joint
        // position only
        mtx44fa_result[2] = k_mtx44f_id;
        mtx44fa_result[2][3][0] = solution->a;
    }


    { // third joint
        const double gamma_complement = solution->gamma - M_PI;
        const double cos_gamma_complement = std::cos(gamma_complement);
        const double sin_gamma_complement = std::sin(gamma_complement);

        MFloatMatrix mtx_local(k_mtx44f_id);

        mtx_local[0][0] = cos_gamma_complement; mtx_local[0][1] = sin_gamma_complement;
        mtx_local[1][0] = -sin_gamma_complement; mtx_local[1][1] = cos_gamma_complement;

        mtx44fa_result[3] = mtx_local;
    }

    { // fourth joint
        MFloatMatrix mtx_local = k_mtx44f_id;

        mtx_local[3][0] = solution->b;

        mtx44fa_result[4] = mtx_local;
    }

    { // fifth joint
        MFloatMatrix mtx_angle_comp(k_mtx44f_id); // realigns the last joint to the basis
        mtx_angle_comp[0][0] = +solution->cos_alpha; mtx_angle_comp[0][1] = +solution->sin_alpha;
        mtx_angle_comp[1][0] = -solution->sin_alpha; mtx_angle_comp[1][1] = +solution->cos_alpha;
        // Transpose here is a cheap inversion because we only care
        //   about the rotation part of the matrix and we know it to be orthonormal
        mtx44fa_result[5] = mtx_handle * mtx_basis.transpose() * mtx_angle_comp;

        mtx44fa_result[5][3][0] = 0.f;
        mtx44fa_result[5][3][1] = 0.f;
        mtx44fa_result[5][3][2] = 0.f;
    }
}


/// \brief
///
/// \param solution a pointer to the triangle data containing the solution to convert to matrices
/// \param mtx_basis Matrix containing the basis (atV and upV transform) to transform to
/// \param mtx_handle Matrix containing world space of the handle to re-align the last transform to
/// \param mtx44fa_result A sequence of matrices representing transforms
///                         at key points in the triangle and algined to the sides.<br>
///                       This is an output argument.
/// \return This is a state-modifying (side-effect only) function
static
void result_from_triangle_flat(const Triangle* const solution,
                               const MFloatMatrix& mtx_basis,
                               const MFloatMatrix& mtx_handle,
                               std::vector<MFloatMatrix>& mtx44fa_result) // output argument
{
    const MFloatMatrix mtx_inverse_basis = mtx_basis.inverse();
    { // zeroth joint under the root, is ID, will get based in the displacement loop
      // this is going to be the plain IK basis
        mtx44fa_result[0] = mtx_basis;
    }

    { // first joint
        MFloatMatrix mtx_local(k_mtx44f_id);
        mtx_local[0][0] = solution->cos_beta; mtx_local[0][1] = solution->sin_beta;
        mtx_local[1][0] = -solution->sin_beta; mtx_local[1][1] = solution->cos_beta;

        mtx44fa_result[1] = mtx_local * mtx_basis;
    }

    { // second joint
        MFloatMatrix mtx_local(k_mtx44f_id);
        // orientation
        mtx_local[0][0] = solution->cos_beta; mtx_local[0][1] = solution->sin_beta;
        mtx_local[1][0] = -solution->sin_beta; mtx_local[1][1] = solution->cos_beta;
        // position
        mtx_local[3][0] = solution->cos_beta * solution->a;
        mtx_local[3][1] = solution->sin_beta * solution->a;

        mtx44fa_result[2] = mtx_local * mtx_basis;
    }

    // Not scoped super-locally because it's used for a couple joints in a row
    const double gamma_complement = solution->gamma + solution->beta - M_PI;
    const double cos_gamma_complement = std::cos(gamma_complement);
    const double sin_gamma_complement = std::sin(gamma_complement);

    { // third joint
        MFloatMatrix mtx_local(k_mtx44f_id);

        mtx_local[0][0] = cos_gamma_complement; mtx_local[0][1] = sin_gamma_complement;
        mtx_local[1][0] = -sin_gamma_complement; mtx_local[1][1] = cos_gamma_complement;

        // position, same as previous joint
        mtx_local[3][0] = solution->cos_beta * solution->a;
        mtx_local[3][1] = solution->sin_beta * solution->a;

        mtx44fa_result[3] = mtx_local * mtx_basis;
    }

    { // fourth joint
        MFloatMatrix mtx_local = k_mtx44f_id;
        // orientation, same as previous joint
        mtx_local[0][0] = cos_gamma_complement; mtx_local[0][1] = sin_gamma_complement;
        mtx_local[1][0] = -sin_gamma_complement; mtx_local[1][1] = cos_gamma_complement;

        mtx_local[3][0] = solution->c;

        mtx44fa_result[4] = mtx_local * mtx_basis;
    }

    { // fifth joint
        mtx44fa_result[5] = mtx_handle;

        mtx44fa_result[5][3][0] = mtx_basis[3][0] + mtx_basis[0][0] * solution->c;
        mtx44fa_result[5][3][1] = mtx_basis[3][1] + mtx_basis[0][1] * solution->c;
        mtx44fa_result[5][3][2] = mtx_basis[3][2] + mtx_basis[0][2] * solution->c;
    }
}


MStatus VChain_soft::compute(const MPlug & plug, MDataBlock & data_block)
{
    if (plug != out_result)
    {
        return MStatus::kUnknownParameter;
    }

    // ----- Cartesian Basis (Aim + UpV) -----
    const MDataHandle hdl_root = data_block.inputValue(in_root);
    const MDataHandle hdl_handle = data_block.inputValue(in_handle);
    const MDataHandle hdl_up_vector = data_block.inputValue(in_up_vector);

    const MFloatMatrix mtx_root = hdl_root.asFloatMatrix();
    const MFloatMatrix mtx_handle = hdl_handle.asFloatMatrix();
    const MFloatMatrix mtx_up_vector = hdl_up_vector.asFloatMatrix();

    const MFloatVector pos_root(mtx_root[3]);
    const MFloatVector pos_handle(mtx_handle[3]);
    const MFloatVector pos_up_vector(mtx_up_vector[3]);

    //    X
    const MFloatVector x_direction = pos_handle - pos_root;
    const float x_distance = x_direction.length();
    const MFloatVector x_cartesian = x_direction.normal();
    //    Y
    const MFloatVector upv_direction = pos_up_vector - pos_handle;
    const MFloatVector y_direction = ((x_cartesian * upv_direction) * x_cartesian) - upv_direction;
    const MFloatVector y_cartesian = y_direction.normal();
    //    Z
    const MFloatVector z_cartesian = x_cartesian ^ y_cartesian;

    // Basis in matrix form
    const MDataHandle hdl_parent_frame = data_block.inputValue(in_parent_frame);
    const MFloatMatrix mtx_parent_frame = hdl_parent_frame.asFloatMatrix();

    auto mtx_ik_basis = matrix_from_vectors(x_cartesian,
                                            y_cartesian,
                                            z_cartesian,
                                            pos_root) * mtx_parent_frame.inverse();

    // ----- Rigid Triangle -----
    //    length parameters
    const MDataHandle hdl_length1 = data_block.inputValue(in_length1);
    const MDataHandle hdl_length2 = data_block.inputValue(in_length2);

    const MDataHandle hdl_max_compression = data_block.inputValue(in_max_compression);
    const MDataHandle hdl_min_extension = data_block.inputValue(in_min_extension);
    const MDataHandle hdl_max_extension = data_block.inputValue(in_max_extension);

    const float length1 = hdl_length1.asFloat();
    const float length2 = hdl_length2.asFloat();

    const float max_compression = hdl_max_compression.asFloat();
    const float min_extension = hdl_min_extension.asFloat();
    const float max_extension = hdl_max_extension.asFloat();

    const float max_length_by_bones = length1 + length2;

    MDataHandle hdl_solver_choice = data_block.inputValue(in_solver_choice);
    const short solver_choice = hdl_solver_choice.asShort();

    //    length solution
    rigid_tri.a = length1;
    rigid_tri.b = length2;

    // Length dampening doesn't require a rigid solution to be softened, it simply
    //   limits the effector when it's reaching for the handle, part of the rigid solution
    const float d_chain = max_length_by_bones;
    const float d_a = d_chain * min_extension;
    if (solver_choice == k_solver_length_dampening && d_a < x_distance)
    {
        // As found here: http://www.softimageblog.com/archives/108
        // Credit to Andy Nicholas
        // Graphed here https://www.desmos.com/calculator/l7n6nszofy
        const float d_soft = d_chain - d_a;
        const float x = x_distance;

        rigid_tri.c = d_soft * (1.0 - std::pow(M_E, (d_a - x) / d_soft)) + d_a;
    } // end of length dampening
    else
    {
        rigid_tri.c = std::max(std::min(x_distance, max_length_by_bones * min_extension),
                                        max_length_by_bones * max_compression );
    }

    //    angular solution
    //      beta
    rigid_tri.cos_beta = rigid_tri.a < 0.00001f ? 0.f : // defending against a div by 0
        ((rigid_tri.a * rigid_tri.a) + (rigid_tri.c * rigid_tri.c) - (rigid_tri.b * rigid_tri.b)) /
                            (2.f * rigid_tri.a * rigid_tri.c);
    rigid_tri.beta = std::acosf(rigid_tri.cos_beta);
    rigid_tri.sin_beta = std::sin(rigid_tri.beta);

    //      gamma
    rigid_tri.cos_gamma = rigid_tri.c < 0.00001f ? 1.f : // defending against a div by 0
        ((rigid_tri.a * rigid_tri.a) + (rigid_tri.b * rigid_tri.b) - (rigid_tri.c * rigid_tri.c)) /
                            (2.f * rigid_tri.a * rigid_tri.b);
    rigid_tri.gamma = std::acosf(rigid_tri.cos_gamma);
    rigid_tri.sin_gamma = std::sin(rigid_tri.gamma);

    //      alpha
    rigid_tri.alpha = M_PI - rigid_tri.gamma - rigid_tri.beta;
    rigid_tri.sin_alpha = std::sin(rigid_tri.alpha);
    rigid_tri.cos_alpha = std::cos(rigid_tri.alpha);

    // ----- Soft Triangle -----
    //  Start with the rigid solution, we'll switch if the extension demands it
    Triangle* chosen_solution = &rigid_tri;

    // We only care for the soft solution if this is true.
    if (rigid_tri.c < x_distance)
    {
        switch (solver_choice)
        {
            case k_solver_rigid :
            {
                // take no action...
            } break;

            case k_solver_height_reduction :
            {
                //  - length solution
                const float max_extended_length = max_length_by_bones * max_extension;
                soft_tri.c = std::min(x_distance, max_extended_length);

                const float stretch_allowance = max_extended_length - rigid_tri.c;
                // Div by 0 should not be an issue given the values taking us to this branch
                // Also this ratio allows trivial addition of non-linear response
                const float current_stretch_ratio = (soft_tri.c - rigid_tri.c) / stretch_allowance;

                const float height_of_rigid_tri = rigid_tri.sin_beta * rigid_tri.a;
                const float height_of_soft_tri = height_of_rigid_tri * (1.f - current_stretch_ratio);

                const float c1 = ((rigid_tri.cos_beta * rigid_tri.a) / rigid_tri.c) * soft_tri.c;
                const float c2 = ((rigid_tri.cos_alpha * rigid_tri.b) / rigid_tri.c) * soft_tri.c;

                soft_tri.a = std::sqrtf((c1 * c1) + (height_of_soft_tri * height_of_soft_tri));
                soft_tri.b = std::sqrtf((c2 * c2) + (height_of_soft_tri * height_of_soft_tri));
                // @maybe: b is probably unused and therefore an optimization point

                //  - angular solution
                soft_tri.beta = std::atanf(height_of_soft_tri / c1);
                soft_tri.sin_beta = std::sin(soft_tri.beta);
                soft_tri.cos_beta = std::cos(soft_tri.beta);

                soft_tri.alpha = std::atanf(height_of_soft_tri / c2);
                soft_tri.sin_alpha = std::sin(soft_tri.alpha);
                soft_tri.cos_alpha = std::cos(soft_tri.alpha);

                soft_tri.gamma = M_PI - soft_tri.alpha - soft_tri.beta;
                soft_tri.sin_gamma = std::sin(soft_tri.gamma);
                soft_tri.cos_gamma = std::cos(soft_tri.gamma);

                chosen_solution = &soft_tri;
            } break;

            case k_solver_height_lock :
            {
                //  - length solution
                const float max_extended_length = max_length_by_bones * max_extension;
                soft_tri.c = std::min(x_distance, max_extended_length);

                const float stretch_allowance = max_extended_length - rigid_tri.c;

                const float height_of_rigid_tri = rigid_tri.sin_beta * rigid_tri.a;

                const float c1 = ((rigid_tri.cos_beta * rigid_tri.a) / rigid_tri.c) * soft_tri.c;
                const float c2 = ((rigid_tri.cos_alpha * rigid_tri.b) / rigid_tri.c) * soft_tri.c;

                soft_tri.a = std::sqrtf((c1 * c1) + (height_of_rigid_tri * height_of_rigid_tri));
                soft_tri.b = std::sqrtf((c2 * c2) + (height_of_rigid_tri * height_of_rigid_tri));
                // @maybe: b is probably unused and therefore an optimization point

                //  - angular solution
                soft_tri.beta = std::atanf(height_of_rigid_tri / c1);
                soft_tri.sin_beta = std::sin(soft_tri.beta);
                soft_tri.cos_beta = std::cos(soft_tri.beta);

                soft_tri.alpha = std::atanf(height_of_rigid_tri / c2);
                soft_tri.sin_alpha = std::sin(soft_tri.alpha);
                soft_tri.cos_alpha = std::cos(soft_tri.alpha);

                soft_tri.gamma = M_PI - soft_tri.alpha - soft_tri.beta;
                soft_tri.sin_gamma = std::sin(soft_tri.gamma);
                soft_tri.cos_gamma = std::cos(soft_tri.gamma);

                chosen_solution = &soft_tri;
            } break;

            case k_solver_length_dampening :
            {
                // Already addressed in rigid solution.
                // Technically this case shouldn't occur since if dampening is on
                //   the opening conditional over-scoping this would be false, but precision
                //   issues might still get us here.
            } break;

            default:
            {
                // There isn't some magic default and we shouldn't get here.
                // It's still going to work as the rigid solution is already configured.
                // @maybe: We might want to error in this case.
                //         Reaching default would mean the enum was pushed out of range somehow.
            }
        }
    }

    //  ------- Result -------
    // Transpose here is a cheap inversion because we only care
    //   about the rotation part of the matrix and we know it to be orthonormal
    const bool hrc_mode = data_block.inputValue(in_hierarchical_mode).asBool();

    if (hrc_mode)
    {
        result_from_triangle_hierarchical(chosen_solution,
                                          mtx_ik_basis,
                                          mtx_handle,
                                          mtx44fa_result); // output argument
    }
    else
    {
        result_from_triangle_flat(chosen_solution,
                                  mtx_ik_basis,
                                  mtx_handle,
                                  mtx44fa_result); // output argument
    }


    //  ------- Output -------
    MArrayDataHandle hdl_result = data_block.outputArrayValue(out_result);
    uint32_t element_count = hdl_result.elementCount();
    uint32_t i = 0;

    // Loop over all elements in the output array that also exist in the
    //   local storage. This could leave you with less than the full
    //   array covered if you don't have all outputs connected, but that's not
    //   really an issue. Most of the time that happens while connections are being made
    for (; i < element_count && i < mtx44fa_result.size(); ++i)
    {
        hdl_result.jumpToArrayElement(i);

        MDataHandle hdl_result_element = hdl_result.outputValue();

        hdl_result_element.setMFloatMatrix(mtx44fa_result[i]);
    }

    // Nothing prevents the user from connecting more plugs than the solver will
    //   actually populate, this will roll over those and set them to id.
    for (; i < element_count; ++i)
    {
        hdl_result.jumpToArrayElement(i);

        MDataHandle hdl_result_element = hdl_result.outputValue();

        hdl_result_element.setMFloatMatrix(k_mtx44f_id);
    }

    hdl_result.setAllClean();

    return MStatus::kSuccess;
}

