#include <iostream>
#include <vector>
#include <utility>

#ifndef BF_H
#define BF_H

bool sort_pair_compare(std::pair<double, int> a, std::pair<double, int> b)
{
  return a.first < b.first;
}

class fixCenterBF
{
public:
  int num;
  std::vector<double> distance_vector;
  std::vector<double> fcbf_value;
  
  fixCenterBF(std::vector<double> &distance_vector, int num)
  {
    this->num = num;
    this->distance_vector = std::vector<double>(num);
    for (unsigned int i = 0; i < num; i++)
    {
      this->distance_vector[i] = distance_vector[i];
    }
    this->fcbf_value = std::vector<double>(num);
    this->train();
  }

  void train()
  {
    int num = this->num;
    std::vector<std::pair<double, int> > dist_index_vec(num);
    for (unsigned int i = 0; i < num; i++)
    {
      dist_index_vec[i] = std::make_pair(this->distance_vector[i], i);
    }
    std::sort(dist_index_vec.begin(), dist_index_vec.end(), sort_pair_compare);
    std::vector<int> index(num);
    for (unsigned int i = 0; i < num; i++)
    {
      this->distance_vector[i] = dist_index_vec[i].first;
      index[i] = dist_index_vec[i].second;
    }

    int lastpos = num - 1;
    double lastval = -1.0;
    for (int j = num - 1; j >= 0; j--)
    {
      if (lastval != this->distance_vector[j])
      {
        lastpos = j;
      }
      lastval = this->distance_vector[j];
      this->fcbf_value[index[j]] = ((double)(lastpos + 1)) / ((double)num);
    }
  }

  void get_fitted(std::vector<double> &fit_value)
  {
    for (int i = 0; i < this->num; i++)
    {
      fit_value[i] = this->fcbf_value[i];
    }
  }

  void predict(std::vector<double> &rank_predict, std::vector<double> new_distance, int predict_num)
  {
    std::vector<double> new_distance_vec(predict_num);
    for (int j = 0; j < predict_num; j++)
    {
      new_distance_vec[j] = new_distance[j];
    }
    std::sort(new_distance_vec.begin(), new_distance_vec.end());

    int train_cursor = 0, predict_cursor = 0;
    while (predict_cursor < predict_num && train_cursor < this->num)
    {
      if (new_distance_vec[predict_cursor] > this->distance_vector[train_cursor])
      {
        rank_predict[train_cursor] = (double)predict_cursor;
        train_cursor++;
      }
      else
      {
        predict_cursor++;
      }
    }
    while (train_cursor < this->num)
    {
      rank_predict[train_cursor] = (double)(predict_num);
      train_cursor++;
    }

    for (int j = 0; j < this->num; j++)
    {
      rank_predict[j] = (rank_predict[j]) / predict_num;
    }
  }

  void free_fixCenterBF()
  {
    // free(this->distance_vector);
    // free(this->fcbf_value);
    // this->distance_vector = NULL;
    // this->fcbf_value = NULL;
  }
};

class BF
{

public:
  int fix_center_num;
  int num;
  std::vector<fixCenterBF *> fcbf_vector;

  BF(std::vector<std::vector<double> > &distance_matrix, int fix_center_num, int num)
  {
    this->fix_center_num = fix_center_num;
    this->num = num;
    this->fcbf_vector = std::vector<fixCenterBF *>(fix_center_num);
    for (unsigned int i = 0; i < this->fix_center_num; i++)
    {
      this->fcbf_vector[i] = new fixCenterBF(distance_matrix[i], num);
    }
    this->train();
  }

  void train()
  {
    for (unsigned int i = 0; i < this->fix_center_num; i++)
    {
      (*this->fcbf_vector[i]).train();
    }
  }

  void get_fitted(std::vector<std::vector<double> > &predict_matrix)
  {
    for (unsigned int i = 0; i < this->fix_center_num; i++)
    {
      (*this->fcbf_vector[i]).get_fitted(predict_matrix[i]);
    }
    return;
  }

  /*
   * @ new_distance: a double type matrix with fix_center_num row, 
   * and each row's center matches the center used for trainning.  
   */
  void predict(std::vector<std::vector<double> > &predict_matrix, std::vector<std::vector<double> > &new_distance, int pred_num)
  {
    for (unsigned int i = 0; i < this->fix_center_num; i++)
    {
      (*this->fcbf_vector[i]).predict(predict_matrix[i], new_distance[i], pred_num);
    }
    return;
  }

  void free_BF()
  {
    // for (unsigned int i = 0; i < this->fix_center_num; i++)
    // {
    //   (*this->fcbf_vector[i]).free_fixCenterBF();
    // }
    for (unsigned int i = 0; i < this->fix_center_num; i++)
    {
      delete this->fcbf_vector[i];
    }
    // free(this->fcbf_vector);
    // this->fcbf_vector = NULL;
  }
};

#endif