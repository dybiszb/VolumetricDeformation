#ifndef VOLUMETRICDEFORMATION_TAPERING_H
#define VOLUMETRICDEFORMATION_TAPERING_H

void linearTaperZ(Eigen::MatrixXd& V, Eigen::MatrixXd& N, double a0, double a1) {
    for (int i = 0; i < V.rows(); ++i) {
        double r = a0 + a1 * V(i, 2);
        double rPrim = a1;
        double x = V(i, 0);
        double y = V(i, 1);

        V(i, 0) = r * x;
        V(i, 1) = r* y;

        double n_x = N(i, 0);
        double n_y = N(i, 1);
        double n_z = N(i, 2);

        N(i, 0) = n_x * r;
        N(i, 1) = n_y * r;
        N(i, 2) = -r * rPrim * n_x * x - r * rPrim * n_y * y + r * r * n_z;
    }
}

void inverseLinearTaperZ(Eigen::MatrixXd& V, Eigen::MatrixXd& N, double a0, double a1) {
    for (int i = 0; i < V.rows(); ++i) {
        double r = a0 + a1 * V(i, 2);
        double rPrim = a1;
        double X = V(i, 0);
        double Y = V(i, 1);

        V(i, 0) = X / r;
        V(i, 1) = Y / r;

        // TODO: Calculate Normals
    }
}

#endif
