/*
 * ARC: Additive uncorrelated Relaxed Clock models for BEAST 2.
 *
 * Copyright (C) 2019-2026 Xavier Didelot, Igor Siveroni and Erik Volz.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
package arc.evolution.branchratemodel;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;

/**
 * Tests for the clock.rate (C1) semantics of {@link ARClockModel}.
 *
 * Note that rateProbs is always supplied at full dimension. If it is given at
 * a lower dimension, initAndValidate expands it with Randomizer.nextDouble(),
 * which is correct behaviour for starting an MCMC chain but makes results
 * incomparable between two models built in the same run.
 */
public class ARClockModelTest {

    /** Three taxa, so five nodes and four branches. */
    private static final String NEWICK = "((A:1.0,B:1.0):1.0,C:2.0):0.0;";

    /** Four probabilities, one per branch, so no random initialisation. */
    private static final String PROBS = "0.1 0.3 0.6 0.9";

    private static TreeParser tree() {
        TreeParser t = new TreeParser();
        t.initByName("newick", NEWICK, "IsLabelledNewick", true, "adjustTipHeights", false);
        return t;
    }

    private static RealParameter fixed(String value) {
        RealParameter p = new RealParameter();
        p.initByName("value", value, "estimate", false);
        return p;
    }

    private static ARClockModel model(Object... extraInputs) {
        Object[] base = {
            "tree", tree(),
            "rateProbs", new RealParameter(PROBS),
            "rateMean", new RealParameter("0.005"),
            "rateOmega", new RealParameter("0.5"),
            "numSites", 1000
        };
        Object[] all = new Object[base.length + extraInputs.length];
        System.arraycopy(base, 0, all, 0, base.length);
        System.arraycopy(extraInputs, 0, all, base.length, extraInputs.length);

        ARClockModel m = new ARClockModel();
        m.initByName(all);
        return m;
    }

    private static double[] branchRates(ARClockModel m) {
        TreeParser t = (TreeParser) m.treeInput.get();
        double[] rates = new double[t.getNodeCount() - 1];
        int k = 0;
        for (int i = 0; i < t.getNodeCount(); i++) {
            Node n = t.getNode(i);
            if (!n.isRoot()) {
                rates[k++] = m.getRateForBranch(n);
            }
        }
        return rates;
    }

    /** Omitting clock.rate must work, and must give strictly positive rates. */
    @Test
    public void clockRateIsOptional() {
        double[] rates = branchRates(model());
        assertEquals(4, rates.length);
        for (double r : rates) {
            assertTrue("branch rate should be positive, was " + r, r > 0.0);
        }
    }

    /**
     * clock.rate must scale every branch rate. Before this was fixed the input
     * was read and then discarded, so setting it changed nothing.
     */
    @Test
    public void clockRateScalesEveryBranchRate() {
        double[] plain  = branchRates(model());
        double[] scaled = branchRates(model("clock.rate", fixed("2.0")));

        assertEquals(plain.length, scaled.length);
        for (int i = 0; i < plain.length; i++) {
            assertEquals("branch " + i, 2.0 * plain[i], scaled[i], 1e-12 * scaled[i]);
        }
    }

    /**
     * clock.rate and rateMean are both mean clock rates, so estimating both is
     * exactly non-identifiable and must be refused rather than sampled.
     */
    @Test
    public void estimatingBothMeanRatesIsRefused() {
        try {
            model("clock.rate", new RealParameter("2.0"));
            fail("expected initAndValidate to reject two estimated mean rates");
        } catch (RuntimeException expected) {
            assertTrue("unhelpful message: " + expected.getMessage(),
                       expected.getMessage().contains("non-identifiable"));
        }
    }

    /** Fixing clock.rate leaves the model identifiable, so it must be accepted. */
    @Test
    public void fixedClockRateWithEstimatedRateMeanIsAccepted() {
        double[] rates = branchRates(model("clock.rate", fixed("2.0")));
        for (double r : rates) {
            assertTrue(r > 0.0);
        }
    }
}
