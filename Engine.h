class Engine {
public:
    virtual void apply(Particle& particle, double deltaTime) = 0;
};