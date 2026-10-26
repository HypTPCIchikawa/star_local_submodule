import copy
import unittest

from prepare_tag_study import select_response_variant


class ResponseVariantTests(unittest.TestCase):
    def config(self):
        return dict(numeric=dict(responseAlphaMin=.25, responseAlphaMax=6.),
                    responseVariants=dict(left=dict(model='gaussian_left_exp', alphaStart=1.)))

    def test_explicit_variant(self):
        cfg = self.config()
        select_response_variant(cfg, 'left')
        self.assertEqual(cfg['numeric']['responseModel'], 'gaussian_left_exp')
        self.assertEqual(cfg['numeric']['responseAlphaStart'], 1.)
        self.assertEqual(cfg['selectedResponseVariant'], 'left')

    def test_no_implicit_response(self):
        with self.assertRaises(ValueError):
            select_response_variant(self.config(), None)

    def test_unknown_variant(self):
        with self.assertRaises(ValueError):
            select_response_variant(self.config(), 'right')

    def test_wrong_family(self):
        cfg = self.config()
        cfg['responseVariants']['left']['model'] = 'free_foreign_peak'
        with self.assertRaises(ValueError):
            select_response_variant(cfg, 'left')

    def test_alpha_support(self):
        for alpha in (.25, 6., -1., float('nan'), True):
            cfg = self.config()
            cfg['responseVariants']['left']['alphaStart'] = alpha
            with self.subTest(alpha=alpha), self.assertRaises(ValueError):
                select_response_variant(cfg, 'left')

    def test_no_silent_shape_fields(self):
        cfg = self.config()
        cfg['responseVariants']['left']['foreignMean'] = 22.
        with self.assertRaises(ValueError):
            select_response_variant(cfg, 'left')

    def test_original_gaussian_settings_unchanged(self):
        cfg = dict(numeric=dict(sigmaLogStep=.02))
        expected = copy.deepcopy(cfg)
        select_response_variant(cfg, None)
        self.assertEqual(cfg, expected)

    def test_no_variant_without_config(self):
        with self.assertRaises(ValueError):
            select_response_variant({}, 'left')


if __name__ == '__main__':
    unittest.main()
