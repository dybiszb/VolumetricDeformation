#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>

#include "twisting.h"
#include "scaling.h"
#include "tapering.h"
#include "bending.h"

//==============================================================================
// Loaded Model
//==============================================================================
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd N;

//==============================================================================
// Scaling Factors
//==============================================================================
double scalingFactorX = 1.1;
double scalingFactorY = 1.1;
double scalingFactorZ = 1.1;

//==============================================================================
// Tapering factors
//==============================================================================
double taperingAlpha0 = 1.0;
double taperingAlpha1 = 0.3;

//==============================================================================
// Twisting Factors
//==============================================================================
double twistScale = 90;

//==============================================================================
// Bending Factors
//==============================================================================
double centerOfBend = 0.0;
double yMax = 0.1;
double yMin = -0.1;
double bendRate = 0.2;

//==============================================================================
// General Settings
//==============================================================================
Eigen::Matrix3d axisStart;
Eigen::Matrix3d axisEnd;
double axisScale = 0.1;
double normalsScale = 0.01;
bool showNormals = false;
igl::viewer::Viewer viewer;

void updateScene(igl::viewer::Viewer &viewer) {
    viewer.data.clear();
    viewer.data.set_mesh(V, F);
    viewer.data.set_normals(N);
    if(showNormals) viewer.data.add_edges(V, V + N * normalsScale, N);
    viewer.data.add_edges(axisStart, axisScale * axisEnd, axisEnd);
}


void initializeCoordinatesFrame() {
    axisStart << 0, 0, 0, 0, 0, 0, 0, 0, 0;
    axisEnd << 1, 0, 0, 0, 1, 0, 0, 0, 1;
}

void moveModelToCenterOfCoordinates() {
    Eigen::MatrixXd centerOfMass = V.colwise().sum();
    centerOfMass /= V.rows();
    for (int i = 0; i < V.rows(); ++i) V.row(i) -= centerOfMass;
}

void loadModel() {
    igl::readOFF("./resources/bunny.off", V, F);
    moveModelToCenterOfCoordinates();
    igl::per_vertex_normals(V, F, N);
}

void initializeUI(igl::viewer::Viewer &viewer) {
    // Extend viewer menu
    viewer.callback_init = [&](igl::viewer::Viewer &viewer) {
        viewer.ngui->addWindow(Eigen::Vector2i(1030, -30), "Transformations");

        /* ------------------- */
        /* ----- Scaling ----- */
        /* ------------------- */
        viewer.ngui->addGroup("Scaling");
        viewer.ngui->addVariable("Scale X axis", scalingFactorX);
        viewer.ngui->addVariable("Scale Y axis", scalingFactorY);
        viewer.ngui->addVariable("Scale Z axis", scalingFactorZ);
        viewer.ngui->addButton("Scale", [&]() {
            scaleModel(V, N, scalingFactorX, scalingFactorY, scalingFactorZ);
            updateScene(viewer);
        });

        viewer.ngui->addButton("Inverse Scale", [&]() {
            scaleModel(V, N, 1.0 / scalingFactorX, 1.0 / scalingFactorY,
                       1.0 / scalingFactorZ);
            updateScene(viewer);
        });

        /* --------------------- */
        /* ----- Tapering ----- */
        /* --------------------- */
        viewer.ngui->addGroup("Tapering");
        viewer.ngui->addVariable("Tapering Function [alpha0]",
                                 taperingAlpha0);
        viewer.ngui->addVariable("Tapering Function [alpha1]",
                                 taperingAlpha1);

        viewer.ngui->addButton("Taper", [&]() {
            linearTaperZ(V, N, taperingAlpha0, taperingAlpha1);
            updateScene(viewer);
        });

        viewer.ngui->addButton("Inverse Taper", [&]() {
            inverseLinearTaperZ(V, N, taperingAlpha0, taperingAlpha1);
            updateScene(viewer);
        });

        /* ----------------------- */
        /* ----- Axial Twist ----- */
        /* ----------------------- */
        viewer.ngui->addGroup("Twisting");
        viewer.ngui->addVariable<double>("Strength along Z", [&](double degree) {
            twistScale = degree;
        }, [&]() {
            return twistScale;
        });
        viewer.ngui->addButton("Twist", [&]() {
            twistModelZ(V, N, twistScale);
            updateScene(viewer);
        });
        viewer.ngui->addButton("Inverse Twist", [&]() {
            inverseTwistModelZ(V, N, twistScale);
            updateScene(viewer);
        });

        /* ------------------- */
        /* ----- Bending ------ */
        /* ------------------- */
        viewer.ngui->addGroup("Bending");
        viewer.ngui->addVariable("Center of bend", centerOfBend);
        viewer.ngui->addVariable("Y min", yMin);
        viewer.ngui->addVariable("Y max", yMax);
        viewer.ngui->addVariable("Bend rate", bendRate);

        viewer.ngui->addButton("Bend", [&]() {
            bendModel(V, N, centerOfBend, yMax, yMin, bendRate);
            updateScene(viewer);
        });

        /* ------------------- */
        /* ----- Global ------ */
        /* ------------------- */
        viewer.ngui->addGroup("Global Settings");
        viewer.ngui->addVariable<bool>("Normals Scaling", [&](bool show) {
            showNormals = show;
            updateScene(viewer);
        }, [&]() {
            return showNormals;
        });
        viewer.ngui->addVariable<double>("Normals Scaling", [&](double scale) {
            normalsScale = scale;
            updateScene(viewer);
        }, [&]() {
            return normalsScale;
        });
        viewer.ngui->addButton("Reset Model", [&]() {
            loadModel();
            updateScene(viewer);
        });

        viewer.screen->performLayout();
        return false;
    };
}

int main(int argc, char *argv[]) {
    loadModel();
    initializeCoordinatesFrame();
    initializeUI(viewer);
    updateScene(viewer);
    viewer.launch();
}




////////////////////// LEGACY

//double taperingAlpha2 = 0.3;

//enum TaperingAxis {
//    X = 0, Y, Z,
//} taperingAxis = Y;
//enum TaperingType {
//    Linear = 0, Quadratic
//} taperingType = Linear;


//        viewer.ngui->addVariable("Tapering Function [alpha2]",
//                                 taperingAlpha2);
//        viewer.ngui->addVariable<TaperingType>("Function Type",
//                                               taperingType)->setItems(
//                {"Linear", "Quadratic"});
//        viewer.ngui->addVariable<TaperingAxis>("Fixed Axis",
//                                               taperingAxis)->setItems(
//                {"X", "Y", "Z"});
//


//
//using TaperingStrategyT = std::function<void(double, double &, double &)>;
//
//void linearTaper(double fixedAxisValue, double &r, double &rPrim) {
//    r = taperingAlpha0 + taperingAlpha1 * fixedAxisValue;
//    rPrim = taperingAlpha1;
//}
//
//void inverseLinearTaper(double fixedAxisValue, double &r, double &rPrim) {
//    r = 1.0 / (taperingAlpha0 + taperingAlpha1 * fixedAxisValue);
////    rPrim = -taperingAlpha1 * (1.0 / (fixedAxisValue * fixedAxisValue));
////    rPrim = - 2 * (1.0 / (r * r));
//}
//
//void quadraticTaper(double fixedAxisValue, double &r, double &rPrim) {
//    r = taperingAlpha0 + taperingAlpha1 * fixedAxisValue +
//        taperingAlpha2 * fixedAxisValue * fixedAxisValue;
//    rPrim = taperingAlpha1 + 2 * taperingAlpha2 * fixedAxisValue;
//}
//
//void inverseQuadraticTaper(double fixedAxisValue, double &r, double &rPrim) {
//    r = 1.0 / (taperingAlpha0 + taperingAlpha1 * fixedAxisValue +
//               taperingAlpha2 * fixedAxisValue * fixedAxisValue);
////    rPrim = -taperingAlpha1 * (1.0 / (fixedAxisValue * fixedAxisValue)) -
////            2.0 * taperingAlpha2 *
////            (1.0 / (fixedAxisValue * fixedAxisValue * fixedAxisValue));
//
////    rPrim = - 2 * (1.0 / (r * r));
//}
//
//void taperModel(const TaperingStrategyT &taperingF) {
//    double fixedAxis, coord1, coord2,fixedAxisCoord, r, rPrim;
//    for (int i = 0; i < V.rows(); ++i) {
//
//        // Get Appropriate Coordinates
//        switch (taperingAxis) {
//            case X:
//                fixedAxisCoord = 0;
//                fixedAxis = V(i, fixedAxisCoord);
//                coord1 = 1;
//                coord2 = 2;
//                break;
//            case Y:
//                fixedAxisCoord = 1;
//                fixedAxis = V(i, fixedAxisCoord);
//                coord1 = 0;
//                coord2 = 2;
//                break;
//            case Z:
//                fixedAxisCoord = 2;
//                fixedAxis = V(i, fixedAxisCoord);
//                coord1 = 0;
//                coord2 = 1;
//                break;
//            default:
//                throw std::runtime_error(
//                        "Wrong axis specified during tapering");
//        }
//
//        taperingF(fixedAxis, r, rPrim);
//
//        // Update Vertices
//        V(i, coord1) *= r;
//        V(i, coord2) *= r;
//
//        // Update Normals
//        // TODO:
////        N(i, coord1) *=r;
////        N(i, coord2) *=r;
////        N(i, fixedAxisCoord) = -r * rPrim * coord1 * N(i, coord1)
////                -r * rPrim * coord2 * N(i, coord2) + r * r * N(i, fixedAxisCoord);
//    }
//}
