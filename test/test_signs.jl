using Random: Xoshiro

@testset "Sign patterns" begin
    s = SignPattern([-1 -1; 1 0]; labels = [:prey, :predator])
    @test s.labels == ["prey", "predator"]
    @test quasi_sign_stability(s) isa NamedTuple
    @test quasi_sign_stability(Xoshiro(1), s; samples = 2000).fraction == 1      # sign-stable
    p = press_sign_determinacy(Xoshiro(1), s; samples = 2000)
    @test p.stable == 1
    @test p.sign == [0 -1; 1 1]          # pressing the predator lowers prey; pressing prey raises the predator
    @test p.sign isa Matrix{Int8} && all(p.determinacy .== 1)
    @test p.labels == ["prey", "predator"]
    @test Set(applicable_metrics(s)) == Set([quasi_sign_stability, press_sign_determinacy])
    @test_throws ArgumentError SignPattern([2 0; 0 -1])
    @test_throws DimensionMismatch SignPattern([-1 0; 0 -1]; labels = ["a"])
    @test_throws DimensionMismatch SignPattern([-1 0; 0 -1]; unknown = falses(3, 3))

    # Without unknown links, a draw consumes the random stream exactly as S .* rand(...).
    a, b = Xoshiro(7), Xoshiro(7)
    @test ResilienceCore._draw(a, s) == s.S .* rand(b, 2, 2)
    @test rand(a) == rand(b)

    # A link of unknown sign takes both signs, with the drawn magnitude.
    u = SignPattern([-1 0; 0 -1]; unknown = Bool[0 1; 0 0])
    @test u.S[1, 2] == 0
    draws = [ResilienceCore._draw(Xoshiro(k), u)[1, 2] for k in 1:400]
    @test all(x -> 0 < abs(x) < 1, draws)
    @test 0.4 < count(>(0), draws) / 400 < 0.6
    # J = [-a b; ±c -d]: stable always when the unknown link is negative (det = ad + bc), and
    # when positive only if ad > bc (probability 1/2), so 3/4 overall
    v = SignPattern([-1 1; 0 -1]; unknown = Bool[0 0; 1 0])
    q = quasi_sign_stability(Xoshiro(3), v; samples = 4000).fraction
    @test abs(q - 0.75) < 0.03
    w = SignPattern([0 1; 0 0]; unknown = Bool[0 0; 1 0])
    qw = quasi_sign_stability(Xoshiro(3), w; samples = 4000).fraction
    @test qw == 0                         # zero diagonal: never strictly stable
    m = SignPattern([-1 1; 0 -1]; unknown = Bool[0 0; 1 0])
    z = SignPattern([-1 1; 1 -1])
    @test quasi_sign_stability(Xoshiro(5), z; samples = 4000).fraction <
          quasi_sign_stability(Xoshiro(5), m; samples = 4000).fraction
end
