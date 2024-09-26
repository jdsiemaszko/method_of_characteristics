from fluidPoint import GenericFlowElement

class BoundaryLine(GenericFlowElement):
    def __init__(self, p, q, v_plus, v_minus, name="Generic Boundary"):
        """
        defines a boundary line segment with pre-defined invariants
        p: point
        q: point
        """
        self.name = name
        self.endpoints = (p, q)
        super.__init__(v_plus, v_minus)



