# Brute-force check of: greedy with l steps covers >= (1-(1-1/k)^l) * OPT_k  (max coverage)
import itertools, random, math
random.seed(1)
worst = 1e9
for trial in range(3000):
    n = random.randint(6, 14); m = random.randint(4, 9); k = random.randint(1, 3)
    sets = [frozenset(random.sample(range(n), random.randint(1, n//2+1))) for _ in range(m)]
    opt = max(len(frozenset().union(*c)) for c in itertools.combinations(sets, k))
    for l in range(1, m+1):
        cov = set(); avail = list(range(m))
        for _ in range(l):
            g = [(len(sets[i]-cov), -i) for i in avail]       # ties: earliest index
            best = max(g)
            if best[0] <= 0: break
            i = -best[1]; avail.remove(i); cov |= sets[i]
        bound = (1-(1-1/k)**l) * opt
        worst = min(worst, len(cov) - bound)
        assert len(cov) >= bound - 1e-9, (sets, k, l, len(cov), opt)
        assert 1-(1-1/k)**l >= 1-math.exp(-l/k) - 1e-12
print("all instances satisfy the bound; min slack", round(worst, 4))
