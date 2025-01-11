#include "stdafx.h"
#include "Include/Utils/rotation3.h"
#include "platform_exception.h"
#include "Solver.h"

namespace platform { namespace math {
	Solver::Solver(const std::map<EMotorPosition, MotorInfo>& motors)
		: joint_(0,0,0)
		, motors_(motors){	
	}

	bool Solver::CalculateAngle(const PlatformState& state, std::map<EMotorPosition, double>& angles){
		PlatformState state1;
		state1.heihgt = state.heihgt;
		state1.roll = state.roll;
		state1.tangage = state.tangage;

		calculateCrankPoints(state1);
		bool tmp = true;
		angles[LEFT] = calculateMotorAngle(motors_[LEFT], left_crank_pl_, state1);
		angles[RIGHT] = calculateMotorAngle(motors_[RIGHT], right_crank_pl_, state1);
		angles[CENTER] = calculateMotorAngle(motors_[CENTER], center_crank_pl_, state1);
		/*fprintf( stdout, "LEFT:%.3f; RIGHT:%.3f; CENTER:%.3f\n", angles[LEFT], angles[RIGHT], angles[CENTER]);
		fflush( stdout);*/
		return tmp;
	}

	void Solver::calculateCrankPoints(const PlatformState& state) {
		
		la::rotation_3 rot_tangage(state.tangage, la::point_3(0,1,0));
		la::rotation_3 rot_roll(state.roll, la::point_3(1,0,0));

		joint_.z = state.heihgt;

		left_crank_pl_ = rot_roll * (rot_tangage * motors_[LEFT].connect_pos) + joint_;
		right_crank_pl_ = rot_roll * (rot_tangage * motors_[RIGHT].connect_pos) + joint_;
		center_crank_pl_ = rot_roll * (rot_tangage * motors_[CENTER].connect_pos) + joint_;
	}

	bool Solver::calculateMotorAngle(const MotorInfo& info, const la::point_3& platform_con, const PlatformState& state, double& out)
	{
		double x = la::abs(info.shaft_pos - platform_con);
		double A = x * x + info.connecting_rod_lenght * info.connecting_rod_lenght - info.crank_lenght * info.crank_lenght;
	//		double A = std::pow(la::abs(info.shaft_pos - platform_con), 2) + std::pow(info.connecting_rod_lenght, 2) - std::pow(info.crank_lenght, 2);
		double B = -2 * info.connecting_rod_lenght * (platform_con.y - info.shaft_pos.y);
		double C = 2 * info.connecting_rod_lenght * (platform_con.z - info.shaft_pos.z);
				
		double discr = B * B + C * C - A * A;
		if(discr < 0)
		{
			fprintf( stdout, "discr < 0\n");
			fflush( stdout);
			return false;
		}
		
		double angle = -2*atan((B + sqrt(discr))/(A - C));
	
		out = angle;
		return true;
	}
	
}}
