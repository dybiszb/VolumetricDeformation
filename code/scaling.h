#ifndef VOLUMETRICDEFORMATION_SCALING_H
#define VOLUMETRICDEFORMATION_SCALING_H

void scaleModel(Eigen::MatrixXd& V, Eigen::MatrixXd& N, double a1, double a2, double a3) {
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

#endif
