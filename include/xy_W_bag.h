#ifndef xy_W_bag_H
#define xy_W_bag_H

#include <vector>

class xy_W_bag
{
public:
    xy_W_bag();

    void W_in(const std::vector<std::vector<double>> &data);
    const std::vector<std::vector<double>> Get_W() const {return bag_W;};

    void b_in(const std::vector<double> &data);
    const std::vector<double> Get_b() const {return bag_b;};

    void d_W_in(const std::vector<std::vector<double>> &data);
    const std::vector<std::vector<double>> Get_d_W() const {return bag_d_W;};

    void d_b_in(const std::vector<double> &data);
    const std::vector<double> Get_d_b() const {return bag_d_b;};

    int Wb_comb();
    const std::vector<std::vector<double>> Get_Wb() const {return bag_Wb;};

private:
    std::vector<std::vector<double>> bag_W;
    std::vector<double> bag_b;
    std::vector<std::vector<double>> bag_d_W;
    std::vector<double> bag_d_b;
    std::vector<std::vector<double>> bag_Wb;


};
#endif
