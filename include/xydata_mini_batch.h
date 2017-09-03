#ifndef xydata_mini_batch_H
#define xydata_mini_batch_H

#include <vector>
#include <string>
#include "xy_W_bag.h"

class xydata_mini_batch
{
public:
    xydata_mini_batch();

    void data_in(const std::vector<std::vector<double>> &data);
    const std::vector<std::vector<double>> Get_data() const {return data_batch;};

    void targ_in(const std::vector<std::vector<double>> &data);
    const std::vector<std::vector<double>> Get_targ() const {return targ_batch;};

    void data_shuffle();

    xy_W_bag Get_mini_batch(unsigned int &id);
    xy_W_bag Get_mini_batch_t(unsigned int &id);

    void make_mini_batch(unsigned int &nb);

    unsigned int Get_nbatch() const {return mini_batch.size();};

private:
    std::vector<std::vector<double>> data_batch;
    std::vector<xy_W_bag> mini_batch;

    std::vector<std::vector<double>> targ_batch;
    std::vector<xy_W_bag> mini_batch_t;

};
#endif
