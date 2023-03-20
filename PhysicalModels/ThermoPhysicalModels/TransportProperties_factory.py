
class FluidTransportProperty:
    def __init__(self) -> None:
        pass

    def getDynamicViscosity(self, sol):
        raise NotImplementedError()

    def getThermalConductivity(self, sol):
        raise NotImplementedError()

class ConstantFluidTransportProperty(FluidTransportProperty):
    def __init__(self, mu: float, kappa: float):
        super().__init__()
        self.mu = mu
        self.kappa = kappa

    def getDynamicViscosity(self, sol):
        del sol
        return self.mu

    def getThermalConductivity(self, sol):
        del sol
        return self.kappa 
    