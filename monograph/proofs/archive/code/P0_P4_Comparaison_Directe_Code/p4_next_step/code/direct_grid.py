"""Certified finite-state baseline; state is one grid index, without a tree.

Uses the compiled Model and its geometric certificate from the P4 audit.
All arithmetic is rational/integer, with size bounded independently of time.
The output is an exact rational Bernoulli parameter (no float conversion).
"""
from fractions import Fraction as F
from certificate import ceil, ceil_log2, point


class DirectGrid:
    def __init__(self, model):
        self.model = model
        self.matrices = tuple(tuple(F(v, model.H) for v in m) for m in model.mats)
        a, b, c, d = self.matrices[0]
        self.alpha, self.beta = a + c - b - d, b + d
        assert self.alpha != 0
        vals = [self.alpha * x + self.beta for x in model.cert.J]
        self.lo, self.hi = min(vals), max(vals)
        self.epsilon = F(1, 1 << model.p)
        self.N = max(1, ceil((self.hi - self.lo) / ((1-model.cert.rho)*self.epsilon)))
        self.spacing = (self.hi-self.lo)/self.N
        self.B = ceil_log2(self.N+1)
        self.error_grid = self.spacing / (2*(1-model.cert.rho))
        self.error_model = model.cert.parameter_TV_constant * model.parameter_error_bound
        self.error_total = self.error_grid + self.error_model
        assert self.error_grid <= self.epsilon/2
        assert self.error_model <= self.epsilon/8
        assert self.error_total <= 5*self.epsilon/8
        self.initial = self.quantize(self.alpha*F(1, 2)+self.beta)

    def output(self, state):
        if type(state) is not int or not 0 <= state <= self.N:
            raise ValueError('invalid grid state')
        return self.lo + state*self.spacing

    def quantize(self, q):
        z = (q-self.lo)/self.spacing
        # Ties go to the larger index. Clamping cannot increase error on Q.
        i = (2*z.numerator+z.denominator)//(2*z.denominator)
        return min(self.N, max(0, i))

    def continuous_update(self, q, symbol):
        if type(symbol) is not int or symbol not in (0, 1):
            raise ValueError('non-binary symbol')
        x = (q-self.beta)/self.alpha
        return self.alpha*point(self.matrices[symbol], x)+self.beta

    def advance(self, state, symbol):
        return self.quantize(self.continuous_update(self.output(state), symbol))


class HybridCodec:
    """Choose once at compilation; the choice is static, not history-dependent."""
    def __init__(self, model):
        from codec import Codec
        self.model = model
        grid = DirectGrid(model)
        fractal = Codec(model)
        self.grid_bits, self.fractal_bits = grid.B, fractal.B
        self.mode = 'grid' if grid.B <= fractal.B else 'fractal'
        self.engine = grid if self.mode == 'grid' else fractal
        self.B = self.engine.B
        self.initial = grid.initial if self.mode == 'grid' else 0

    def advance(self, state, symbol):
        return self.engine.advance(state, symbol)

    def output(self, state):
        if self.mode == 'grid':
            return self.engine.output(state)
        return F(self.engine.output(state), self.model.Q)
