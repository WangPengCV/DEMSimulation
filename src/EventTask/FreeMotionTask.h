#include "EventTask.h"

class FreeMotionTask : public EventTask {
public:
    FreeMotionTask(Eigen::Vector3d& position, Eigen::Vector3d& velocity);
    virtual ~FreeMotionTask() override = default;

    const Eigen::Vector3d&  getPosition() const; 
    const Eigen::Vector3d&  getVelocity() const; 


private:
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;

};