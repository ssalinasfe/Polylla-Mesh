double get_min_angle_non_zero(int *poly, int len_poly, double *r);
double get_angle_three_points(double *r, int i, int j, int k);
double get_max_edge(int *poly, int len_poly, double *r);
double get_min_edge(int *poly, int len_poly, double *r);
//double get_signed_area_poly(int *poly, int length_poly, double *r);
void get_metrics_of_polygon(double *metrics, int *poly, int len_poly, double *r);

//int GreaterEqualthan(double a, double b, double epsilon);
//int Equality(double a, double b, double epsilon);
//double dist(double x0, double y0, double x1, double y1);