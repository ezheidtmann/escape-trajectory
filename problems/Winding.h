#ifndef R3WINDING_H
#define R3WINDING_H

#include <cvode.h>
#include <math.h>

/**
 * Utility class to count winding number for the RestrictedThree problem.
 *
 * The winding number is the number of radians through which the projectile
 * has traveled around a certain center of rotation: If the projectile 
 * takes a loopy path, this number should be large. If it makes a straight
 * shot out of the system or performs a number of switchbacks, this number
 * should be small.
 */
class Winding {
  public:

	inline void addStep(realtype x, realtype y);

	/**
	 * Get the current winding number
	 */
	realtype getWindingNumber() { return _windingNumber; }
		
	Winding():
		_windingNumber(0),
		_last_defined(false)
	{}

  private:
	realtype _windingNumber;
	realtype _last_x, _last_y;
	bool _last_defined;
};

/**
 * Add a step of motion. If the projectile moves from the point 
 * \f$\overrightarrow{r_0} = \langle x_0, y_0 \rangle\f$ to 
 * \f$\overrightarrow{r} = \langle x, y \rangle\f$, then
 * the magnitude of the angle between the vectors is found by
 * 
 * \f{eqnarray*}
 * 	\left|\Theta\right| &=& \arccos \left( \hat{r_0} \cdot \hat{r} \right)\\
 * 	&=& \arccos \left( \frac{xx_0 + yy_0}
 * 		{\sqrt{\left(x^2+y^2\right)\left(x_0^2+y_0^2\right)}} \right)
 * \f}
 *
 * The direction of rotation is given by the direction of the cross product
 * of the vectors. Because both \f$\overrightarrow{r_0}\f$ and
 * \f$\overrightarrow{r}\f$ are 2-dimensional, the cross product is
 * zero in all components except the \f$z\f$ component, which has value
 * \f$x_0y - xy_0\f$.
 *
 * If this quantity is greater than zero, \f$\Theta\f$ is positive 
 * (corresponding to a counter-clockwise rotation). If
 * the quantity is less than zero, \f$\Theta\f$ is negative (corresponding
 * to a clockwise rotation).
 *
 * If the computed angle is invalid, as might occur near the origin or when
 * some of the values are very small, we define \f$\Theta \equiv 0\f$.
 *
 * When called for the first time, only the previous point is updated.
 * No computation is performed.
 *
 * @param x X-coordinate, with origin at center of rotation
 * @param y Y-coordinate, with origin at center of rotation
 */
void Winding::addStep(realtype x, realtype y)
{
	static realtype theta;

	if (_last_defined)
	{
		theta = acos((x*_last_x + y*_last_y) / 
			sqrt( (x*x + y*y) * (_last_x*_last_x + _last_y*_last_y) )
		);
		theta *= (_last_x*y - x*_last_y) > 0 ? 1 : -1;
	
		// Only increment winding number if theta is sane.
		if (theta == theta) // IEEE equivalent to theta != NaN
			_windingNumber += theta;
		
	}
	
	_last_x = x;
	_last_y = y;
	_last_defined = true;
}

#endif

