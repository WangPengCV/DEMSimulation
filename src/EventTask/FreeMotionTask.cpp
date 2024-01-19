#include "FreeMotionTask.h"

FreeMotionTask::FreeMotionTask(Eigen::Vector3d& position, Eigen::Vector3d& velocity)
    :position(position),velocity(velocity)
{

}

const Eigen::Vector3d &FreeMotionTask::getPosition() const
{
    return position;
}
const Eigen::Vector3d &FreeMotionTask::getVelocity() const
{
    return velocity;
}