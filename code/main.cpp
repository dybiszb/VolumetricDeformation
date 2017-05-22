#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd N;

double scalingFactorX = 1.1;
double scalingFactorY = 1.1;
double scalingFactorZ = 1.1;

enum TaperingAxis {
    X = 0, Y, Z,
} taperingAxis = Y;
enum TaperingType {
    Linear = 0, Quadratic
} taperingType = Linear;

double taperingAlpha0 = 1.0;
double taperingAlpha1 = 0.3;
double taperingAlpha2 = 0.3;

double twistScale = 0.1;

Eigen::Matrix3d axisStart;
Eigen::Matrix3d axisEnd;
double axisScale = 0.1;
double normalsScale = 0.05;
bool showNormals = false;

igl::viewer::Viewer viewer;

void updateScene(igl::viewer::Viewer &viewer) {
    viewer.data.clear();
    viewer.data.set_mesh(V, F);
    viewer.data.set_normals(N);
    if(showNormals) viewer.data.add_edges(V, V + N * normalsScale, N);
    viewer.data.add_edges(axisStart, axisScale * axisEnd, axisEnd);
}

void scaleModel(double a1, double a2, double a3) {
    // Vertices
    Eigen::MatrixXd m(3, 3);
    m << a1, 0, 0, 0, a2, 0, 0, 0, a3;
    V = V * m.transpose();

    // Normals
    Eigen::MatrixXd invTJac(3, 3);
    invTJac << 1.0 / a1, 0, 0, 0, 1.0 / a2, 0, 0, 0, 1.0 / a3;
    N = (invTJac * N.transpose()).transpose();
    N.rowwise().normalize();
}

using TaperingStrategyT = std::function<void(double, double &, double &)>;

void linearTaper(double fixedAxisValue, double &r, double &rPrim) {
    r = taperingAlpha0 + taperingAlpha1 * fixedAxisValue;
    rPrim = taperingAlpha1;
}

void inverseLinearTaper(double fixedAxisValue, double &r, double &rPrim) {
    r = 1.0 / (taperingAlpha0 + taperingAlpha1 * fixedAxisValue);
//    rPrim = -taperingAlpha1 * (1.0 / (fixedAxisValue * fixedAxisValue));
//    rPrim = - 2 * (1.0 / (r * r));
}

void quadraticTaper(double fixedAxisValue, double &r, double &rPrim) {
    r = taperingAlpha0 + taperingAlpha1 * fixedAxisValue +
        taperingAlpha2 * fixedAxisValue * fixedAxisValue;
    rPrim = taperingAlpha1 + 2 * taperingAlpha2 * fixedAxisValue;
}

void inverseQuadraticTaper(double fixedAxisValue, double &r, double &rPrim) {
    r = 1.0 / (taperingAlpha0 + taperingAlpha1 * fixedAxisValue +
               taperingAlpha2 * fixedAxisValue * fixedAxisValue);
//    rPrim = -taperingAlpha1 * (1.0 / (fixedAxisValue * fixedAxisValue)) -
//            2.0 * taperingAlpha2 *
//            (1.0 / (fixedAxisValue * fixedAxisValue * fixedAxisValue));

//    rPrim = - 2 * (1.0 / (r * r));
}

void taperModel(const TaperingStrategyT &taperingF) {
    double fixedAxis, coord1, coord2,fixedAxisCoord, r, rPrim;
    for (int i = 0; i < V.rows(); ++i) {

        // Get Appropriate Coordinates
        switch (taperingAxis) {
            case X:
                fixedAxisCoord = 0;
                fixedAxis = V(i, fixedAxisCoord);
                coord1 = 1;
                coord2 = 2;
                break;
            case Y:
                fixedAxisCoord = 1;
                fixedAxis = V(i, fixedAxisCoord);
                coord1 = 0;
                coord2 = 2;
                break;
            case Z:
                fixedAxisCoord = 2;
                fixedAxis = V(i, fixedAxisCoord);
                coord1 = 0;
                coord2 = 1;
                break;
            default:
                throw std::runtime_error(
                        "Wrong axis specified during tapering");
        }

        taperingF(fixedAxis, r, rPrim);

        // Update Vertices
        V(i, coord1) *= r;
        V(i, coord2) *= r;

        // Update Normals
        // TODO:
//        N(i, coord1) *=r;
//        N(i, coord2) *=r;
//        N(i, fixedAxisCoord) = -r * rPrim * coord1 * N(i, coord1)
//                -r * rPrim * coord2 * N(i, coord2) + r * r * N(i, fixedAxisCoord);
    }
}

void twistModel() {
    for (int i = 0; i < V.rows(); ++i) {
        double z = V(i, 2);
        double theta = 0.0174532925 *  1000 * z;
        double C_theta = cos(theta);
        double S_theta = sin(theta);

        V(i, 0) = V(i, 0) * C_theta - V(i, 1) * S_theta;
        V(i, 1) = V(i, 0) * C_theta + V(i, 1) * S_theta;
    }
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

void initializeUI(igl::viewer::Viewer &viewer) {
    // Extend viewer menu
    viewer.callback_init = [&](igl::viewer::Viewer &viewer) {
        viewer.ngui->addWindow(Eigen::Vector2i(1000, 10), "Transformations");

        /* ------------------- */
        /* ----- Scaling ----- */
        /* ------------------- */
        viewer.ngui->addGroup("Scaling");
        viewer.ngui->addVariable("Scale X axis", scalingFactorX);
        viewer.ngui->addVariable("Scale Y axis", scalingFactorY);
        viewer.ngui->addVariable("Scale Z axis", scalingFactorZ);
        viewer.ngui->addButton("Scale", [&]() {
            scaleModel(scalingFactorX, scalingFactorY, scalingFactorZ);
            updateScene(viewer);
        });

        viewer.ngui->addButton("Inverse Scale", [&]() {
            scaleModel(1.0 / scalingFactorX, 1.0 / scalingFactorY,
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
        viewer.ngui->addVariable("Tapering Function [alpha2]",
                                 taperingAlpha2);
        viewer.ngui->addVariable<TaperingType>("Function Type",
                                               taperingType)->setItems(
                {"Linear", "Quadratic"});
        viewer.ngui->addVariable<TaperingAxis>("Fixed Axis",
                                               taperingAxis)->setItems(
                {"X", "Y", "Z"});
        viewer.ngui->addButton("Taper", [&]() {
            if (taperingType == TaperingType::Linear) {
                taperModel(linearTaper);
            } else {
                taperModel(quadraticTaper);
            }

            updateScene(viewer);
        });

        viewer.ngui->addButton("Inverse Taper", [&]() {
            if (taperingType == TaperingType::Linear) {
                taperModel(inverseLinearTaper);
            } else {
                taperModel(inverseQuadraticTaper);
            }

            updateScene(viewer);
        });

        /* ----------------------- */
        /* ----- Axial Twist ----- */
        /* ----------------------- */
        viewer.ngui->addGroup("Twisting");
        viewer.ngui->addVariable<double>("Degrees (temp)", [&](double degree) {
            twistScale = 0.0174532925 * degree;
        }, [&]() {
            return twistScale;
        });
        viewer.ngui->addButton("Twist", [&]() {
            twistModel();
            updateScene(viewer);
        });




        /* ------------------- */
        /* ----- Normals ----- */
        /* ------------------- */
        viewer.ngui->addGroup("Normals");
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

        viewer.screen->performLayout();
        return false;
    };
}

int main(int argc, char *argv[]) {
    igl::readOFF("./resources/lion.off", V, F);
    moveModelToCenterOfCoordinates();
    igl::per_vertex_normals(V, F, N);

    initializeCoordinatesFrame();
    initializeUI(viewer);
    updateScene(viewer);
    viewer.launch();
}



