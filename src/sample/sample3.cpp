#include "iostream"
#include "cmath"
#include "../lib/lbfgs.h"

using std::clog;
using std::endl;

class TEST_FUNC
{
public:
	TEST_FUNC();
	~TEST_FUNC();
	static lbfgsfloatval_t _Func(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g,
		const int n, const lbfgsfloatval_t step)
	{
		return reinterpret_cast<TEST_FUNC*>(instance)->Func(x, g, n, step);
	}

	lbfgsfloatval_t Func(const lbfgsfloatval_t *x, lbfgsfloatval_t *g,
		const int n, const lbfgsfloatval_t step);

	static int _Progress(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
		const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step,
		int n, int k, int ls)
	{
		return reinterpret_cast<TEST_FUNC*>(instance)->Progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
	}

	int Progress(const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
		const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step,
		int n, int k, int ls);

	int Routine();
private:
	lbfgsfloatval_t *m_x;
};

TEST_FUNC::TEST_FUNC()
{
	m_x = NULL;
	m_x = lbfgs_malloc(3);
	m_x[0] = m_x[1] = m_x[2] = 0.0;
}

TEST_FUNC::~TEST_FUNC()
{
	if (m_x != NULL) lbfgs_free(m_x);
}

// test functions
// 3 = 3*x1 + x2 + 2*x3*x3
// 1 = -3*x1 + 5*x2*x2 + 2*x1*x3
// -12 = 25*x1*x2 + 20*x3
lbfgsfloatval_t TEST_FUNC::Func(const lbfgsfloatval_t *x, lbfgsfloatval_t *g,
		const int n, const lbfgsfloatval_t step)
{
	double f0,f1,f2,temp;
	f0 = 3*x[0] + x[1] + 2*x[2]*x[2] - 3.012; //这里添加一点噪声
	f1 = -3*x[0] + 5*x[1]*x[1] + 2*x[0]*x[2] - 1.0521;
	f2 = 25*x[0]*x[1] + 20*x[2] + 12.10231;
	temp = sqrt(f0*f0+f1*f1+f2*f2);

	g[0] = 0.5*(6*f0+2*f1*(2*x[2]-3)+50*f2*x[1])/temp;
	g[1] = 0.5*(2*f0+20*f1*x[1]+50*f2*x[0])/temp;
	g[2] = 0.5*(8*f0*x[2]+4*f1*x[0]+40*f2)/temp;
	return temp;
}

int TEST_FUNC::Progress(const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
		const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step,
		int n, int k, int ls)
{
	clog << "iteration times: " << k << " fx = " << fx << " gnorm/xnorm = " << gnorm/xnorm << endl;
	clog << x[0] << " " << x[1] << " " << x[2] << endl;

	if (fx < 1e-10) return 1; //这里我们设置一个方程组的整体目标函数值作为终止条件 因为此方程组在0值处梯度不为0 无法使用梯度条件
	return 0;
}

int TEST_FUNC::Routine()
{
	lbfgsfloatval_t fx;

	lbfgs_parameter_t self_para;
	lbfgs_parameter_init(&self_para);
	//self_para.min_step = 1e-30;
	//self_para.max_linesearch = 40;
	//self_para.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;

	int ret = lbfgs(3, m_x, &fx, _Func, _Progress, this, &self_para);

	clog << "L-BFGS optimization terminated with status: " << endl << lbfgs_strerror(ret) << endl;
	clog << m_x[0] << " " << m_x[1] << " " << m_x[2] << endl;
	return ret;
}

int main(int argc, char const *argv[])
{
	TEST_FUNC test;
	test.Routine();
	return 0;
}