#ifndef VOLUMETRICDEFORMATION_TWISTING_H
#define VOLUMETRICDEFORMATION_TWISTING_H

void
twistModelZ(Eigen::MatrixXd &V, Eigen::MatrixXd &N, const double twistScale) {
    const double degreesToRadians = 0.0174532925;
    for (int i = 0; i < V.rows(); ++i) {
        double z = V(i, 2);
        double theta = degreesToRadians * twistScale * z;
        double C_theta = cos(theta);
        double S_theta = sin(theta);

        // Get old values of x nad y
        double x = V(i, 0);
        double y = V(i, 1);

        // Calculate new coordinates
        V(i, 0) = x * C_theta - y * S_theta;
        V(i, 1) = x * S_theta + y * C_theta;

        // Calculate normals
        double n_x = N(i, 0);
        double n_y = N(i, 1);
        double n_z = N(i, 2);
        N(i, 0) = n_x * C_theta - S_theta * n_y;
        N(i, 1) = n_x * S_theta + C_theta * n_y;
        N(i, 2) = n_x * y * theta - n_y * x * theta + n_z;
    }
}

void inverseTwistModelZ(Eigen::MatrixXd &V, Eigen::MatrixXd &N,
                        const double twistScale) {
    const double degreesToRadians = 0.0174532925;
    for (int i = 0; i < V.rows(); ++i) {
        double Z = V(i, 2);
        double theta = degreesToRadians * twistScale * Z;
        double C_theta = cos(theta);
        double S_theta = sin(theta);

        double X = V(i, 0);
        double Y = V(i, 1);

        V(i, 0) = X * C_theta + Y * S_theta;
        V(i, 1) = -X * S_theta + Y * C_theta;

        // TODO : Calculate Normals
    }
}

#endif
