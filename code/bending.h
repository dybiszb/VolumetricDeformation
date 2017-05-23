#ifndef VOLUMETRICDEFORMATION_BENDING_H
#define VOLUMETRICDEFORMATION_BENDING_H


void
bendModel(Eigen::MatrixXd &V, Eigen::MatrixXd &N, double y0,
          double yMax, double yMin, double k) {

    for (int i = 0; i < V.rows(); ++i) {
        // Collect current coordinates
        double x = V(i, 0);
        double y = V(i, 1);
        double z = V(i, 2);

        // Pick y-hat
        double yHat = y;
        if (yHat < yMin) yHat = yMin;
        else if (yHat > yMax) yHat = yMax;

        // Bending angle
        double theta = k * (yHat - y0);
        double C_theta = cos(theta);
        double S_theta = sin(theta);

        // Calculate Y
        if (y >= yMin && y <= yMax) {
            V(i, 1) = -S_theta * (z - (1.0 / k)) + y0;
        } else if (y < yMin) {
            V(i, 1) = -S_theta * (z - (1.0 / k)) + y0 + C_theta * (y - yMin);
        } else if (y > yMax) {
            V(i, 1) = -S_theta * (z - (1.0 / k)) + y0 + C_theta * (y - yMax);
        }

        // Calculate Z
        if (y >= yMin && y <= yMax) {
            V(i, 2) = C_theta * (z - (1.0 / k)) + (1.0 / k);
        } else if (y < yMin) {
            V(i, 2) = C_theta * (z - (1.0 / k)) + (1.0 / k) + S_theta * (y - yMin);
        } else if (y > yMax) {
            V(i, 2) = C_theta * (z - (1.0 / k)) + (1.0 / k) + S_theta * (y - yMax);
        }

        // Calculate k-hat
        double khat = (y == yHat) ? k : 0.0;

        // Calculate normals
        double n_x = N(i, 0);
        double n_y = N(i, 1);
        double n_z = N(i, 2);
        double locRtExp = (1.0 - khat * z);

        N(i, 0) = n_x * locRtExp;
        N(i, 1) = n_y * C_theta - S_theta * locRtExp * n_z;
        N(i, 2) = n_y * S_theta + C_theta * locRtExp * n_z;
    }

}

#endif
