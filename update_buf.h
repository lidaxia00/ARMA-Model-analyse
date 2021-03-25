#ifndef __UPDATE_BUF_H
#define __UPDATE_BUF_H

typedef struct predict_temp
{
	float predict_position;
	float predict_position_cha;
	float predict_positionAR;
	float predict_err_ar;
	float predict_err_arma;
	float *buf_temp;
	float *e_data;
} predict_s;

extern const float sensor_err[];
extern const float sensor_down[];
extern predict_s predict_temp;

void update_bufs(predict_s *p, float d, float e, int buf_size);

#endif