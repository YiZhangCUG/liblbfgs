/*
 *      C library of Limited memory BFGS (L-BFGS).
 *
 * Copyright (c) 1990, Jorge Nocedal
 * Copyright (c) 2007-2010 Naoaki Okazaki
 * All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/* $Id$ */

#ifndef __LBFGS_H__
#define __LBFGS_H__

#ifdef  __cplusplus //c++编译环境中才会定义__cplusplus (plus就是"+"的意思)
//作用是让 C++ 编译器将 extern "C" 声明的代码当作 C 语言代码处理，可以避免 C++ 因符号修饰导致代码不能和C语言库中的符号进行链接的问题。
extern "C" {
#endif/*__cplusplus*/

// 算法库使用的浮点类型定义
// 结合后面的定义可见我们现在使用的其实就是 double 类型
/*
 * The default precision of floating point values is 64bit (double).
 */
#ifndef LBFGS_FLOAT
#define LBFGS_FLOAT     64
#endif/*LBFGS_FLOAT*/

/*
 * Activate optimization routines for IEEE754 floating point values.
 */
#ifndef LBFGS_IEEE_FLOAT
#define LBFGS_IEEE_FLOAT    1
#endif/*LBFGS_IEEE_FLOAT*/

#if     LBFGS_FLOAT == 32
typedef float lbfgsfloatval_t;

#elif   LBFGS_FLOAT == 64
typedef double lbfgsfloatval_t;

#else
#error "libLBFGS supports single (float; LBFGS_FLOAT = 32) or double (double; LBFGS_FLOAT=64) precision only."

#endif

// 枚举类型 lbfgs()返回值 定义了各种返回值对应的别称
// 返回值为负代表错误
/** 
 * \addtogroup liblbfgs_api libLBFGS API
 * @{
 *
 *  The libLBFGS API.
 */

/**
 * Return values of lbfgs().
 * 
 *  Roughly speaking, a negative value indicates an error.
 */
enum {
    /** L-BFGS reaches convergence. */
    LBFGS_SUCCESS = 0,
    LBFGS_CONVERGENCE = 0,
    LBFGS_STOP, //1
    /** The initial variables already minimize the objective function. */
    LBFGS_ALREADY_MINIMIZED, //2

    /** Unknown error. */
    LBFGSERR_UNKNOWNERROR = -1024,
    /** Logic error. */
    LBFGSERR_LOGICERROR, //-1023
    /** Insufficient memory. */
    LBFGSERR_OUTOFMEMORY, //-1022
    /** The minimization process has been canceled. */
    LBFGSERR_CANCELED,
    /** Invalid number of variables specified. */
    LBFGSERR_INVALID_N,
    /** Invalid number of variables (for SSE) specified. */
    LBFGSERR_INVALID_N_SSE,
    /** The array x must be aligned to 16 (for SSE). */
    LBFGSERR_INVALID_X_SSE,
    /** Invalid parameter lbfgs_parameter_t::epsilon specified. */
    LBFGSERR_INVALID_EPSILON,
    /** Invalid parameter lbfgs_parameter_t::past specified. */
    LBFGSERR_INVALID_TESTPERIOD,
    /** Invalid parameter lbfgs_parameter_t::delta specified. */
    LBFGSERR_INVALID_DELTA,
    /** Invalid parameter lbfgs_parameter_t::linesearch specified. */
    LBFGSERR_INVALID_LINESEARCH,
    /** Invalid parameter lbfgs_parameter_t::max_step specified. */
    LBFGSERR_INVALID_MINSTEP,
    /** Invalid parameter lbfgs_parameter_t::max_step specified. */
    LBFGSERR_INVALID_MAXSTEP,
    /** Invalid parameter lbfgs_parameter_t::ftol specified. */
    LBFGSERR_INVALID_FTOL,
    /** Invalid parameter lbfgs_parameter_t::wolfe specified. */
    LBFGSERR_INVALID_WOLFE,
    /** Invalid parameter lbfgs_parameter_t::gtol specified. */
    LBFGSERR_INVALID_GTOL,
    /** Invalid parameter lbfgs_parameter_t::xtol specified. */
    LBFGSERR_INVALID_XTOL,
    /** Invalid parameter lbfgs_parameter_t::max_linesearch specified. */
    LBFGSERR_INVALID_MAXLINESEARCH,
    /** Invalid parameter lbfgs_parameter_t::orthantwise_c specified. */
    LBFGSERR_INVALID_ORTHANTWISE,
    /** Invalid parameter lbfgs_parameter_t::orthantwise_start specified. */
    LBFGSERR_INVALID_ORTHANTWISE_START,
    /** Invalid parameter lbfgs_parameter_t::orthantwise_end specified. */
    LBFGSERR_INVALID_ORTHANTWISE_END,
    /** The line-search step went out of the interval of uncertainty. */
    LBFGSERR_OUTOFINTERVAL,
    /** A logic error occurred; alternatively, the interval of uncertainty
        became too small. */
    LBFGSERR_INCORRECT_TMINMAX,
    /** A rounding error occurred; alternatively, no line-search step
        satisfies the sufficient decrease and curvature conditions. */
    LBFGSERR_ROUNDING_ERROR,
    /** The line-search step became smaller than lbfgs_parameter_t::min_step. */
    LBFGSERR_MINIMUMSTEP,
    /** The line-search step became larger than lbfgs_parameter_t::max_step. */
    LBFGSERR_MAXIMUMSTEP,
    /** The line-search routine reaches the maximum number of evaluations. */
    LBFGSERR_MAXIMUMLINESEARCH,
    /** The algorithm routine reaches the maximum number of iterations. */
    LBFGSERR_MAXIMUMITERATION,
    /** Relative width of the interval of uncertainty is at most
        lbfgs_parameter_t::xtol. */
    LBFGSERR_WIDTHTOOSMALL,
    /** A logic error (negative line-search step) occurred. */
    LBFGSERR_INVALIDPARAMETERS,
    /** The current search direction increases the objective function value. */
    LBFGSERR_INCREASEGRADIENT,
};

// 枚举类型 线性搜索方法
// 0 MoreThuente方法
// 1 Armijo条件方法
// 2 标准Wolfe条件方法
// 3 增强Wolfe条件方法
/**
 * Line search algorithms.
 */
enum {
    /** The default algorithm (MoreThuente method). */
    LBFGS_LINESEARCH_DEFAULT = 0,
    /** MoreThuente method proposd by More and Thuente. */
    LBFGS_LINESEARCH_MORETHUENTE = 0,
    /**
     * Backtracking method with the Armijo condition.
     *  The backtracking method finds the step length such that it satisfies
     *  the sufficient decrease (Armijo) condition,
     *    - f(x + a * d) <= f(x) + lbfgs_parameter_t::ftol * a * g(x)^T d,
     *
     *  where x is the current point, d is the current search direction, and
     *  a is the step length.
     */
    LBFGS_LINESEARCH_BACKTRACKING_ARMIJO = 1,
    /** The backtracking method with the defualt (regular Wolfe) condition. */
    LBFGS_LINESEARCH_BACKTRACKING = 2,
    /**
     * Backtracking method with regular Wolfe condition.
     *  The backtracking method finds the step length such that it satisfies
     *  both the Armijo condition (LBFGS_LINESEARCH_BACKTRACKING_ARMIJO)
     *  and the curvature condition,
     *    - g(x + a * d)^T d >= lbfgs_parameter_t::wolfe * g(x)^T d,
     *
     *  where x is the current point, d is the current search direction, and
     *  a is the step length.
     */
    LBFGS_LINESEARCH_BACKTRACKING_WOLFE = 2,
    /**
     * Backtracking method with strong Wolfe condition.
     *  The backtracking method finds the step length such that it satisfies
     *  both the Armijo condition (LBFGS_LINESEARCH_BACKTRACKING_ARMIJO)
     *  and the following condition,
     *    - |g(x + a * d)^T d| <= lbfgs_parameter_t::wolfe * |g(x)^T d|,
     *
     *  where x is the current point, d is the current search direction, and
     *  a is the step length.
     */
    LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE = 3,
};

// L-BFGS参数类型。参数很多，简要说明如下：
// m L-BFGS算法中储存的前序sk与yk向量个数，这个值控制了算法使用的内存多少，默认值为6（不建议小于3的值），值多大近似精度越高，计算量也越大。
// epsilon 迭代的终止精度，默认值为1e-5
// past 以delta(不同迭代次数的目标函数值)为基础的迭代终止条件数，past代表了以多少迭代次数之前的目标函数值作为delta计算的间隔，默认值为0，
//      即不以delta为迭代终止条件。
// delta (f' - f) / f 不同迭代次数时目标函数之差与当前目标函数值之比，但past不为0时会计算。
// max_iterations 最大迭代次数，为0时表示一直迭代到终止条件被满足或出现其他错误。
// linesearch 线性搜索方式，由此文件前述枚举类型定义。
// max_linesearch 每次迭代中线性搜索的最大次数，默认值为40
// min_step 线性搜索中的最小步长，默认值为1e-20
// max_step 线性搜索中的最大步长，默认值为1e+20
// ftol 线性搜索的精度值，默认值为1e-4，取值范围（0-0.5）。
// wolfe Wolfe线性搜索中的控制参数，默认值为0.9，大于ftol小于1.0
// gtol 线性搜索中的控制参数，默认值为0.9，大于ftol小于1.0
// xtol 浮点数精度，默认值为1e-16
// orthantwise_c 模型参数x的L1模的乘积参数，默认值为0.0，此时算法即为L2模形式，当此参数大于0时，算法即为OWL-QN
// orthantwise_start 开始计算模型参数x的L1模的迭代序号
// orthantwise_end 终止计算模型参数x的L1模的迭代序号
/**
 * L-BFGS optimization parameters.
 *  Call lbfgs_parameter_init() function to initialize parameters to the
 *  default values.
 */
typedef struct {
    /**
     * The number of corrections to approximate the inverse hessian matrix.
     *  The L-BFGS routine stores the computation results of previous \ref m
     *  iterations to approximate the inverse hessian matrix of the current
     *  iteration. This parameter controls the size of the limited memories
     *  (corrections). The default value is \c 6. Values less than \c 3 are
     *  not recommended. Large values will result in excessive computing time.
     */
    int             m;

    /**
     * Epsilon for convergence test.
     *  This parameter determines the accuracy with which the solution is to
     *  be found. A minimization terminates when
     *      ||g|| < \ref epsilon * max(1, ||x||),
     *  where ||.|| denotes the Euclidean (L2) norm. The default value is
     *  \c 1e-5.
     */
    lbfgsfloatval_t epsilon;

    /**
     * Distance for delta-based convergence test.
     *  This parameter determines the distance, in iterations, to compute
     *  the rate of decrease of the objective function. If the value of this
     *  parameter is zero, the library does not perform the delta-based
     *  convergence test. The default value is \c 0.
     */
    int             past;

    /**
     * Delta for convergence test.
     *  This parameter determines the minimum rate of decrease of the
     *  objective function. The library stops iterations when the
     *  following condition is met:
     *      (f' - f) / f < \ref delta,
     *  where f' is the objective value of \ref past iterations ago, and f is
     *  the objective value of the current iteration.
     *  The default value is \c 1e-5.
     */
    lbfgsfloatval_t delta;

    /**
     * The maximum number of iterations.
     *  The lbfgs() function terminates an optimization process with
     *  ::LBFGSERR_MAXIMUMITERATION status code when the iteration count
     *  exceedes this parameter. Setting this parameter to zero continues an
     *  optimization process until a convergence or error. The default value
     *  is \c 0.
     */
    int             max_iterations;

    /**
     * The line search algorithm.
     *  This parameter specifies a line search algorithm to be used by the
     *  L-BFGS routine.
     */
    int             linesearch;

    /**
     * The maximum number of trials for the line search.
     *  This parameter controls the number of function and gradients evaluations
     *  per iteration for the line search routine. The default value is \c 40.
     */
    int             max_linesearch;

    /**
     * The minimum step of the line search routine.
     *  The default value is \c 1e-20. This value need not be modified unless
     *  the exponents are too large for the machine being used, or unless the
     *  problem is extremely badly scaled (in which case the exponents should
     *  be increased).
     */
    lbfgsfloatval_t min_step;

    /**
     * The maximum step of the line search.
     *  The default value is \c 1e+20. This value need not be modified unless
     *  the exponents are too large for the machine being used, or unless the
     *  problem is extremely badly scaled (in which case the exponents should
     *  be increased).
     */
    lbfgsfloatval_t max_step;

    /**
     * A parameter to control the accuracy of the line search routine.
     *  The default value is \c 1e-4. This parameter should be greater
     *  than zero and smaller than \c 0.5.
     */
    lbfgsfloatval_t ftol;

    /**
     * A coefficient for the Wolfe condition.
     *  This parameter is valid only when the backtracking line-search
     *  algorithm is used with the Wolfe condition,
     *  ::LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE or
     *  ::LBFGS_LINESEARCH_BACKTRACKING_WOLFE .
     *  The default value is \c 0.9. This parameter should be greater
     *  the \ref ftol parameter and smaller than \c 1.0.
     */
    lbfgsfloatval_t wolfe;

    /**
     * A parameter to control the accuracy of the line search routine.
     *  The default value is \c 0.9. If the function and gradient
     *  evaluations are inexpensive with respect to the cost of the
     *  iteration (which is sometimes the case when solving very large
     *  problems) it may be advantageous to set this parameter to a small
     *  value. A typical small value is \c 0.1. This parameter shuold be
     *  greater than the \ref ftol parameter (\c 1e-4) and smaller than
     *  \c 1.0.
     */
    lbfgsfloatval_t gtol;

    /**
     * The machine precision for floating-point values.
     *  This parameter must be a positive value set by a client program to
     *  estimate the machine precision. The line search routine will terminate
     *  with the status code (::LBFGSERR_ROUNDING_ERROR) if the relative width
     *  of the interval of uncertainty is less than this parameter.
     */
    lbfgsfloatval_t xtol;

    /**
     * Coeefficient for the L1 norm of variables.
     *  This parameter should be set to zero for standard minimization
     *  problems. Setting this parameter to a positive value activates
     *  Orthant-Wise Limited-memory Quasi-Newton (OWL-QN) method, which
     *  minimizes the objective function F(x) combined with the L1 norm |x|
     *  of the variables, {F(x) + C |x|}. This parameter is the coeefficient
     *  for the |x|, i.e., C. As the L1 norm |x| is not differentiable at
     *  zero, the library modifies function and gradient evaluations from
     *  a client program suitably; a client program thus have only to return
     *  the function value F(x) and gradients G(x) as usual. The default value
     *  is zero.
     */
    lbfgsfloatval_t orthantwise_c;

    /**
     * Start index for computing L1 norm of the variables.
     *  This parameter is valid only for OWL-QN method
     *  (i.e., \ref orthantwise_c != 0). This parameter b (0 <= b < N)
     *  specifies the index number from which the library computes the
     *  L1 norm of the variables x,
     *      |x| := |x_{b}| + |x_{b+1}| + ... + |x_{N}| .
     *  In other words, variables x_1, ..., x_{b-1} are not used for
     *  computing the L1 norm. Setting b (0 < b < N), one can protect
     *  variables, x_1, ..., x_{b-1} (e.g., a bias term of logistic
     *  regression) from being regularized. The default value is zero.
     */
    int             orthantwise_start;

    /**
     * End index for computing L1 norm of the variables.
     *  This parameter is valid only for OWL-QN method
     *  (i.e., \ref orthantwise_c != 0). This parameter e (0 < e <= N)
     *  specifies the index number at which the library stops computing the
     *  L1 norm of the variables x,
     */
    int             orthantwise_end;
} lbfgs_parameter_t;

// 目标函数与其梯度值计算的回调函数模版，参数简要说明如下：
// instance 运行实例的指针，帮助程序正确定位回调函数在内存中的位置。
// x 当前的模型参数值的指针
// g 当前模型参数值对应的梯度指针
// n 模型参数的数量
// step 当前线性搜索所使用的步长
// retval 当前模型参数的目标函数值
/**
 * Callback interface to provide objective function and gradient evaluations.
 *
 *  The lbfgs() function call this function to obtain the values of objective
 *  function and its gradients when needed. A client program must implement
 *  this function to evaluate the values of the objective function and its
 *  gradients, given current values of variables.
 *  
 *  @param  instance    The user data sent for lbfgs() function by the client.
 *  @param  x           The current values of variables.
 *  @param  g           The gradient vector. The callback function must compute
 *                      the gradient values for the current variables.
 *  @param  n           The number of variables.
 *  @param  step        The current step of the line search routine.
 *  @retval lbfgsfloatval_t The value of the objective function for the current
 *                          variables.
 */
typedef lbfgsfloatval_t (*lbfgs_evaluate_t)(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    );

// 进程函数的回调函数模版，参数简要说明如下：
// instance 运行实例的指针，帮助程序正确定位回调函数在内存中的位置。
// x 当前的模型参数值的指针
// g 当前模型参数值对应的梯度指针
// fx 目标函数的值
// xnorm 模型参数数组的L2模长
// gnorm 模型梯度数组的L2模长
// step 当前线性搜索所使用的步长
// n 模型参数的数量
// k 迭代的次数
// ls 此次迭代所使用的线性搜索次数
// retval 返回0则lbfgs()函数继续，否则终止
/**
 * Callback interface to receive the progress of the optimization process.
 *
 *  The lbfgs() function call this function for each iteration. Implementing
 *  this function, a client program can store or display the current progress
 *  of the optimization process.
 *
 *  @param  instance    The user data sent for lbfgs() function by the client.
 *  @param  x           The current values of variables.
 *  @param  g           The current gradient values of variables.
 *  @param  fx          The current value of the objective function.
 *  @param  xnorm       The Euclidean norm of the variables.
 *  @param  gnorm       The Euclidean norm of the gradients.
 *  @param  step        The line-search step used for this iteration.
 *  @param  param       这是我们添加了一个指针以使用参数类型来监控迭代流程
 *  @param  n           The number of variables.
 *  @param  k           The iteration count.
 *  @param  ls          The number of evaluations called for this iteration.
 *  @retval int         Zero to continue the optimization process. Returning a
 *                      non-zero value will cancel the optimization process.
 */
typedef int (*lbfgs_progress_t)(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    const lbfgs_parameter_t param,
    int n,
    int k,
    int ls
    );

/*
A user must implement a function compatible with ::lbfgs_evaluate_t (evaluation
callback) and pass the pointer to the callback function to lbfgs() arguments.
Similarly, a user can implement a function compatible with ::lbfgs_progress_t
(progress callback) to obtain the current progress (e.g., variables, function
value, ||G||, etc) and to cancel the iteration process if necessary.
Implementation of a progress callback is optional: a user can pass \c NULL if
progress notification is not necessary.

// 这里提出了算法使用中的两个要求
// 1. 变量的数量必须是16的倍数
// 2. 变量的储存以16对齐
// 还不太明白为什么要这么要求。这里需要以后再注解。
// 注解：貌似只有使用SEE
In addition, a user must preserve two requirements:
    - The number of variables must be multiples of 16 (this is not 4).
    - The memory block of variable array ::x must be aligned to 16.

This algorithm terminates an optimization
when:

    ||G|| < \epsilon \cdot \max(1, ||x||) .

In this formula, ||.|| denotes the Euclidean norm.
*/

// 下面是L-BFGS的主函数，各个参数的说明简要翻译如下：
// n 数组的长度，也就是待求的模型参数的数量
// x 模型参数数组的指针，函数通过指针直接操作模型数组，所以不需要返回计算结果。一开始赋给函数的数组即为
//   初始模型，函数结束后即为最优化结果
// ptr_fx 目标函数的值的指针，设计成指针可以方便在函数外部监控迭代过程的收敛情况
// proc_evaluate 计算目标函数与目标函数相对于模型的导数的函数名称。这个函数在实际使用中需要用户按照
//               算法库指定的参数形式自行定义。函数功能即为计算非线性最优化问题的目标函数值与相应的
//               模型梯度。L-BFGS将利用这两个量来计算雅各比与近似海森矩阵，进而确定迭代的方向与步长。
// proc_progress 接收迭代过程指标参数的函数名称，这个函数的作为即方便用户以自定义的形式为迭代过程进行
//               监控或显示。
// instance 函数执行时的实例对象，将被proc_evaluate函数与proc_progress函数接收。这个变量的存在
//          是由于如果回调函数是类的成员函数时无法直接调用，需要通过定义一个静态的类成员函数作为新的回调函数
//          来调用回调函数。但是类内的函数指针调用只能通过实例进行，所以必然需要一个空的指针来指向运行时的实例。
//          而这个指针也就一路来到了这里，它的作用就是在运行时帮助程序确定回调函数在内存的位置以保证正确的调用。
// param L-BFGS算法参数类型的指针，指向一个包含算法运行需要的参数的结构体。
// retval 返回值。无错即为0，非0值代表此文件上部枚举类型中的对应错误。此文件下部定义的错误信息显示即利用此返回值与
//        预定义的枚举类型输出相应的错误信息。
/**
 * Start a L-BFGS optimization.
 *
 *  @param  n           The number of variables.
 *  @param  x           The array of variables. A client program can set
 *                      default values for the optimization and receive the
 *                      optimization result through this array. This array
 *                      must be allocated by ::lbfgs_malloc function
 *                      for libLBFGS built with SSE/SSE2 optimization routine
 *                      enabled. The library built without SSE/SSE2
 *                      optimization does not have such a requirement.
 *  @param  ptr_fx      The pointer to the variable that receives the final
 *                      value of the objective function for the variables.
 *                      This argument can be set to \c NULL if the final
 *                      value of the objective function is unnecessary.
 *  @param  proc_evaluate   The callback function to provide function and
 *                          gradient evaluations given a current values of
 *                          variables. A client program must implement a
 *                          callback function compatible with \ref
 *                          lbfgs_evaluate_t and pass the pointer to the
 *                          callback function.
 *  @param  proc_progress   The callback function to receive the progress
 *                          (the number of iterations, the current value of
 *                          the objective function) of the minimization
 *                          process. This argument can be set to \c NULL if
 *                          a progress report is unnecessary.
 *  @param  instance    A user data for the client program. The callback
 *                      functions will receive the value of this argument.
 *  @param  param       The pointer to a structure representing parameters for
 *                      L-BFGS optimization. A client program can set this
 *                      parameter to \c NULL to use the default parameters.
 *                      Call lbfgs_parameter_init() function to fill a
 *                      structure with the default values.
 *  @retval int         The status code. This function returns zero if the
 *                      minimization process terminates without an error. A
 *                      non-zero value indicates an error.
 */
int lbfgs(
    int n,
    lbfgsfloatval_t *x,
    lbfgsfloatval_t *ptr_fx,
    lbfgs_evaluate_t proc_evaluate,
    lbfgs_progress_t proc_progress,
    void *instance,
    lbfgs_parameter_t *param
    );

// 将一个参数类型内的全部变量值重置为默认值
// 如果害怕参数被自己调乱了可以用这个函数将参数重置
/**
 * Initialize L-BFGS parameters to the default values.
 *
 *  Call this function to fill a parameter structure with the default values
 *  and overwrite parameter values if necessary.
 *
 *  @param  param       The pointer to the parameter structure.
 */
void lbfgs_parameter_init(lbfgs_parameter_t *param);

// 开辟浮点类型的数组空间
/**
 * Allocate an array for variables.
 *
 *  This function allocates an array of variables for the convenience of
 *  ::lbfgs function; the function has a requreiemt for a variable array
 *  when libLBFGS is built with SSE/SSE2 optimization routines. A user does
 *  not have to use this function for libLBFGS built without SSE/SSE2
 *  optimization.
 *  
 *  @param  n           The number of variables.
 */
lbfgsfloatval_t* lbfgs_malloc(int n);

// 释放浮点类型的数组空间
/**
 * Free an array of variables.
 *  
 *  @param  x           The array of variables allocated by ::lbfgs_malloc
 *                      function.
 */
void lbfgs_free(lbfgsfloatval_t *x);

// 使用lbfgs()函数的返回值被返回一个包含相应错误信息的字符串。
// 锦上添花，不是必要的部分。
/**
 * Get string description of an lbfgs() return code.
 *
 *  @param err          A value returned by lbfgs().
 */
const char* lbfgs_strerror(int err);

/** @} */
// 这里链接头文件开始位置的声明，保证了整个头文件都会以C的形式被编译
#ifdef  __cplusplus
}
#endif/*__cplusplus*/



/**
@mainpage libLBFGS: a library of Limited-memory Broyden-Fletcher-Goldfarb-Shanno (L-BFGS)

@section intro Introduction

This library is a C port of the implementation of Limited-memory
Broyden-Fletcher-Goldfarb-Shanno (L-BFGS) method written by Jorge Nocedal.
The original FORTRAN source code is available at:
http://www.ece.northwestern.edu/~nocedal/lbfgs.html

The L-BFGS method solves the unconstrainted minimization problem,

<pre>
    minimize F(x), x = (x1, x2, ..., xN),
</pre>

only if the objective function F(x) and its gradient G(x) are computable. The
well-known Newton's method requires computation of the inverse of the hessian
matrix of the objective function. However, the computational cost for the
inverse hessian matrix is expensive especially when the objective function
takes a large number of variables. The L-BFGS method iteratively finds a
minimizer by approximating the inverse hessian matrix by information from last
m iterations. This innovation saves the memory storage and computational time
drastically for large-scaled problems.

Among the various ports of L-BFGS, this library provides several features:
- <b>Optimization with L1-norm (Orthant-Wise Limited-memory Quasi-Newton
  (OWL-QN) method)</b>:
  In addition to standard minimization problems, the library can minimize
  a function F(x) combined with L1-norm |x| of the variables,
  {F(x) + C |x|}, where C is a constant scalar parameter. This feature is
  useful for estimating parameters of sparse log-linear models (e.g.,
  logistic regression and maximum entropy) with L1-regularization (or
  Laplacian prior).
- <b>Clean C code</b>:
  Unlike C codes generated automatically by f2c (Fortran 77 into C converter),
  this port includes changes based on my interpretations, improvements,
  optimizations, and clean-ups so that the ported code would be well-suited
  for a C code. In addition to comments inherited from the original code,
  a number of comments were added through my interpretations.
- <b>Callback interface</b>:
  The library receives function and gradient values via a callback interface.
  The library also notifies the progress of the optimization by invoking a
  callback function. In the original implementation, a user had to set
  function and gradient values every time the function returns for obtaining
  updated values.
- <b>Thread safe</b>:
  The library is thread-safe, which is the secondary gain from the callback
  interface.
- <b>Cross platform.</b> The source code can be compiled on Microsoft Visual
  Studio 2010, GNU C Compiler (gcc), etc.
- <b>Configurable precision</b>: A user can choose single-precision (float)
  or double-precision (double) accuracy by changing ::LBFGS_FLOAT macro.
- <b>SSE/SSE2 optimization</b>:
  This library includes SSE/SSE2 optimization (written in compiler intrinsics)
  for vector arithmetic operations on Intel/AMD processors. The library uses
  SSE for float values and SSE2 for double values. The SSE/SSE2 optimization
  routine is disabled by default.

This library is used by:
- <a href="http://www.chokkan.org/software/crfsuite/">CRFsuite: A fast implementation of Conditional Random Fields (CRFs)</a>
- <a href="http://www.chokkan.org/software/classias/">Classias: A collection of machine-learning algorithms for classification</a>
- <a href="http://www.public.iastate.edu/~gdancik/mlegp/">mlegp: an R package for maximum likelihood estimates for Gaussian processes</a>
- <a href="http://infmath.uibk.ac.at/~matthiasf/imaging2/">imaging2: the imaging2 class library</a>

@section download Download

- <a href="https://github.com/downloads/chokkan/liblbfgs/liblbfgs-1.10.tar.gz">Source code</a>
- <a href="https://github.com/chokkan/liblbfgs">GitHub repository</a>

libLBFGS is distributed under the term of the
<a href="http://opensource.org/licenses/mit-license.php">MIT license</a>.

@section modules Third-party modules
- <a href="http://cran.r-project.org/web/packages/lbfgs/index.html">lbfgs: Limited-memory BFGS Optimization (a wrapper for R)</a> maintained by Antonio Coppola.
- <a href="http://search.cpan.org/~laye/Algorithm-LBFGS-0.16/">Algorithm::LBFGS - Perl extension for L-BFGS</a> maintained by Lei Sun.
- <a href="http://www.cs.kuleuven.be/~bernd/yap-lbfgs/">YAP-LBFGS (an interface to call libLBFGS from YAP Prolog)</a> maintained by Bernd Gutmann.

@section changelog History
- Version 1.10 (2010-12-22):
    - Fixed compiling errors on Mac OS X; this patch was kindly submitted by
      Nic Schraudolph.
    - Reduced compiling warnings on Mac OS X; this patch was kindly submitted
      by Tamas Nepusz.
    - Replaced memalign() with posix_memalign().
    - Updated solution and project files for Microsoft Visual Studio 2010.
- Version 1.9 (2010-01-29):
    - Fixed a mistake in checking the validity of the parameters "ftol" and
      "wolfe"; this was discovered by Kevin S. Van Horn.
- Version 1.8 (2009-07-13):
    - Accepted the patch submitted by Takashi Imamichi;
      the backtracking method now has three criteria for choosing the step
      length:
        - ::LBFGS_LINESEARCH_BACKTRACKING_ARMIJO: sufficient decrease (Armijo)
          condition only
        - ::LBFGS_LINESEARCH_BACKTRACKING_WOLFE: regular Wolfe condition
          (sufficient decrease condition + curvature condition)
        - ::LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE: strong Wolfe condition
    - Updated the documentation to explain the above three criteria.
- Version 1.7 (2009-02-28):
    - Improved OWL-QN routines for stability.
    - Removed the support of OWL-QN method in MoreThuente algorithm because
      it accidentally fails in early stages of iterations for some objectives.
      Because of this change, <b>the OW-LQN method must be used with the
      backtracking algorithm (::LBFGS_LINESEARCH_BACKTRACKING)</b>, or the
      library returns ::LBFGSERR_INVALID_LINESEARCH.
    - Renamed line search algorithms as follows:
        - ::LBFGS_LINESEARCH_BACKTRACKING: regular Wolfe condition.
        - ::LBFGS_LINESEARCH_BACKTRACKING_LOOSE: regular Wolfe condition.
        - ::LBFGS_LINESEARCH_BACKTRACKING_STRONG: strong Wolfe condition.
    - Source code clean-up.
- Version 1.6 (2008-11-02):
    - Improved line-search algorithm with strong Wolfe condition, which was
      contributed by Takashi Imamichi. This routine is now default for
      ::LBFGS_LINESEARCH_BACKTRACKING. The previous line search algorithm
      with regular Wolfe condition is still available as
      ::LBFGS_LINESEARCH_BACKTRACKING_LOOSE.
    - Configurable stop index for L1-norm computation. A member variable
      ::lbfgs_parameter_t::orthantwise_end was added to specify the index
      number at which the library stops computing the L1 norm of the
      variables. This is useful to prevent some variables from being
      regularized by the OW-LQN method.
    - A sample program written in C++ (sample/sample.cpp).
- Version 1.5 (2008-07-10):
    - Configurable starting index for L1-norm computation. A member variable
      ::lbfgs_parameter_t::orthantwise_start was added to specify the index
      number from which the library computes the L1 norm of the variables.
      This is useful to prevent some variables from being regularized by the
      OWL-QN method.
    - Fixed a zero-division error when the initial variables have already
      been a minimizer (reported by Takashi Imamichi). In this case, the
      library returns ::LBFGS_ALREADY_MINIMIZED status code.
    - Defined ::LBFGS_SUCCESS status code as zero; removed unused constants,
      LBFGSFALSE and LBFGSTRUE.
    - Fixed a compile error in an implicit down-cast.
- Version 1.4 (2008-04-25):
    - Configurable line search algorithms. A member variable
      ::lbfgs_parameter_t::linesearch was added to choose either MoreThuente
      method (::LBFGS_LINESEARCH_MORETHUENTE) or backtracking algorithm
      (::LBFGS_LINESEARCH_BACKTRACKING).
    - Fixed a bug: the previous version did not compute psuedo-gradients
      properly in the line search routines for OWL-QN. This bug might quit
      an iteration process too early when the OWL-QN routine was activated
      (0 < ::lbfgs_parameter_t::orthantwise_c).
    - Configure script for POSIX environments.
    - SSE/SSE2 optimizations with GCC.
    - New functions ::lbfgs_malloc and ::lbfgs_free to use SSE/SSE2 routines
      transparently. It is uncessary to use these functions for libLBFGS built
      without SSE/SSE2 routines; you can still use any memory allocators if
      SSE/SSE2 routines are disabled in libLBFGS.
- Version 1.3 (2007-12-16):
    - An API change. An argument was added to lbfgs() function to receive the
      final value of the objective function. This argument can be set to
      \c NULL if the final value is unnecessary.
    - Fixed a null-pointer bug in the sample code (reported by Takashi Imamichi).
    - Added build scripts for Microsoft Visual Studio 2005 and GCC.
    - Added README file.
- Version 1.2 (2007-12-13):
    - Fixed a serious bug in orthant-wise L-BFGS.
      An important variable was used without initialization.
- Version 1.1 (2007-12-01):
    - Implemented orthant-wise L-BFGS.
    - Implemented lbfgs_parameter_init() function.
    - Fixed several bugs.
    - API documentation.
- Version 1.0 (2007-09-20):
    - Initial release.

@section api Documentation

- @ref liblbfgs_api "libLBFGS API"

@section sample Sample code

@include sample.c

@section ack Acknowledgements

The L-BFGS algorithm is described in:
    - Jorge Nocedal.
      Updating Quasi-Newton Matrices with Limited Storage.
      <i>Mathematics of Computation</i>, Vol. 35, No. 151, pp. 773--782, 1980.
    - Dong C. Liu and Jorge Nocedal.
      On the limited memory BFGS method for large scale optimization.
      <i>Mathematical Programming</i> B, Vol. 45, No. 3, pp. 503-528, 1989.

The line search algorithms used in this implementation are described in:
    - John E. Dennis and Robert B. Schnabel.
      <i>Numerical Methods for Unconstrained Optimization and Nonlinear
      Equations</i>, Englewood Cliffs, 1983.
    - Jorge J. More and David J. Thuente.
      Line search algorithm with guaranteed sufficient decrease.
      <i>ACM Transactions on Mathematical Software (TOMS)</i>, Vol. 20, No. 3,
      pp. 286-307, 1994.

This library also implements Orthant-Wise Limited-memory Quasi-Newton (OWL-QN)
method presented in:
    - Galen Andrew and Jianfeng Gao.
      Scalable training of L1-regularized log-linear models.
      In <i>Proceedings of the 24th International Conference on Machine
      Learning (ICML 2007)</i>, pp. 33-40, 2007.

Special thanks go to:
    - Yoshimasa Tsuruoka and Daisuke Okanohara for technical information about
      OWL-QN
    - Takashi Imamichi for the useful enhancements of the backtracking method
    - Kevin S. Van Horn, Nic Schraudolph, and Tamas Nepusz for bug fixes

Finally I would like to thank the original author, Jorge Nocedal, who has been
distributing the effieicnt and explanatory implementation in an open source
licence.

@section reference Reference

- <a href="http://www.ece.northwestern.edu/~nocedal/lbfgs.html">L-BFGS</a> by Jorge Nocedal.
- <a href="http://research.microsoft.com/en-us/downloads/b1eb1016-1738-4bd5-83a9-370c9d498a03/default.aspx">Orthant-Wise Limited-memory Quasi-Newton Optimizer for L1-regularized Objectives</a> by Galen Andrew.
- <a href="http://chasen.org/~taku/software/misc/lbfgs/">C port (via f2c)</a> by Taku Kudo.
- <a href="http://www.alglib.net/optimization/lbfgs.php">C#/C++/Delphi/VisualBasic6 port</a> in ALGLIB.
- <a href="http://cctbx.sourceforge.net/">Computational Crystallography Toolbox</a> includes
  <a href="http://cctbx.sourceforge.net/current_cvs/c_plus_plus/namespacescitbx_1_1lbfgs.html">scitbx::lbfgs</a>.
*/

#endif/*__LBFGS_H__*/
