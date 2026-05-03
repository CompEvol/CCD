package ccd.algorithms.regularisation;

public enum CCDRegularisationStrategy {
    AdditiveOne("Additive-1"),
    AdditiveX("Additive-"),
    PriorOne("Prior-1"),
    PriorScaled("Prior-Scaled");

    private final String name;

    CCDRegularisationStrategy(String name) {
        this.name = name;
    }

    public CCDRegularisationStrategy setValue(double value) {
        if (value <= 0) {
            throw new IllegalArgumentException("Additive/scaling value for regularization/smoothing must be positive.");
        }
        return this;
    }

    @Override
    public String toString() {
        return this.name;
    }
}
