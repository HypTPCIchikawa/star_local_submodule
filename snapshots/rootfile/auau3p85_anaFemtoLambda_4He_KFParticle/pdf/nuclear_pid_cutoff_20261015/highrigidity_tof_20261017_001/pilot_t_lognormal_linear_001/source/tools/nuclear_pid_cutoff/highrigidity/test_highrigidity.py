import unittest
import run_highrigidity as run


class Provider:
    @staticmethod
    def current_parameters(runtime,species,rigidity):
        return (10.,1.,'fixture') if rigidity<=5 else (0.,0.,'outside_domain')


class HighRigidityTests(unittest.TestCase):
    def test_competitors_active_before_global_cap(self):
        self.assertTrue(run.gate_parameters({'cuts':{'maxPOverQ':2.5}},Provider,0,3.)[2])
        self.assertFalse(run.gate_parameters({},Provider,3,5.01)[2])

    def test_native_pilot_geometry(self):
        cfg=dict(nativeWidth=.01,rigidityMin=2.,rigidityMax=2.5,pilotRigidity=[2.005,2.245,2.495])
        self.assertEqual(run.selected_bins(cfg,'pilot'),[201,225,250])
        self.assertEqual(run.selected_bins(cfg,'native'),list(range(201,251)))
        cfg['pilotRigidity']=[2.006]
        with self.assertRaises(ValueError): run.selected_bins(cfg,'pilot')

    def test_actual_and_study_tof_windows_are_distinct(self):
        context={'settings':{'tagDMass2Min':3.268773,'tagDMass2Max':3.693147}}
        current={'cuts':{'m2_selection':False,'highRigidityTof':dict(enabled=True,thresholdPOverQ=2.,species={
            'd':dict(mass2Min=3.268773,mass2Max=3.693147),'he4':dict(mass2Min=3.,mass2Max=4.)})}}
        self.assertEqual(run.require_exact_window('d',dict(sample='TAG_D',exactProductionTofWindow=True),current,context),(3.268773,3.693147))
        with self.assertRaises(ValueError):
            run.require_exact_window('he4',dict(sample='TAG_D',exactProductionTofWindow=True),current,context)
        self.assertEqual(run.require_exact_window('he4',dict(sample='TAG_D',exactProductionTofWindow=False,
            studyTofPolicy='same_saved_D_window'),current,context),(3.268773,3.693147))

    def test_target_reoptimized_not_varied_as_prior(self):
        cfg=dict(allFit=dict(targetShapeMode='free_target_mean_sigma'),shapeVariations=['width_plus'],widthRelativeShift=.1)
        component=dict(mean=10.,sigma=2.)
        self.assertEqual(run.common.varied_response(component,'width_plus',cfg,'d','d'),(10.,2.))
        self.assertEqual(run.common.varied_response(component,'width_plus',cfg,'he3','d'),(10.,2.2))


if __name__=='__main__': unittest.main()
