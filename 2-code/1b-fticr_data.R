


icr_report_positive = read.csv("1-data/fticr/DOC 367-420 Positive Report.csv")
icr_report_negative = read.csv("1-data/fticr/DOC 367-420 Negative Report.csv")
sample_key = read.csv("1-data/sample_key.csv", na = "") %>% dplyr::select(-notes)

meta_negative = make_fticr_meta(icr_report = icr_report_negative)$meta2
meta_positive = make_fticr_meta(icr_report = icr_report_positive)$meta2
data_negative = make_fticr_data(icr_report_negative, sample_key)$data_long_trt
data_positive = make_fticr_data(icr_report_positive, sample_key)$data_long_trt

data_negative_longform = make_fticr_data(icr_report_negative, sample_key)$data_long_blank_corrected
data_positive_longform = make_fticr_data(icr_report_positive, sample_key)$data_long_blank_corrected



meta_combined = 
  meta_negative %>% mutate(mode = "negative") %>% 
  bind_rows(meta_positive %>% mutate(mode = "positive")) %>% 
  group_by(formula) %>% 
  mutate(n = n()) %>% 
  mutate(mode = case_when(n == 2 ~ "both", TRUE ~ mode)) %>% 
  distinct()

###############
###############


gg_vk_domains = 
  gg_vankrev(meta_negative, aes(x = OC, y = HC, color = Class))+
  scale_color_manual(values = PNWColors::pnw_palette("Sunset2"))

gg_vankrev(meta_positive, aes(x = OC, y = HC, color = Class))+
  scale_color_manual(values = PNWColors::pnw_palette("Sunset2"))


gg_vankrev(meta_combined, aes(x = OC, y = HC, color = mode))+
  facet_wrap(~mode)

###############
###############




neg_unique = 
  data_negative %>% 
  group_by(formula, region, site, horizon) %>% 
  mutate(n = n())

neg_unique %>% 
  left_join(meta_combined %>% dplyr::select(formula, HC, OC)) %>% 
  filter(! horizon %in% "B") %>% 
  gg_vankrev(aes(x = OC, y = HC, color = as.character(n)))+
  facet_wrap(~region + n)


neg_unique %>% 
  left_join(meta_combined %>% dplyr::select(formula, HC, OC)) %>% 
  filter(! horizon %in% "B") %>% 
  filter(n == 1) %>% 
  gg_vankrev(aes(x = OC, y = HC, color = transect))+
  facet_wrap(~region + transect)


pos_unique = 
  data_positive %>% 
  group_by(formula, region, site, horizon) %>% 
  mutate(n = n())

pos_unique %>% 
  left_join(meta_combined %>% dplyr::select(formula, HC, OC)) %>% 
  filter(! horizon %in% "B") %>% 
  gg_vankrev(aes(x = OC, y = HC, color = as.character(n)))+
  facet_wrap(~region + n)


pos_unique %>% 
  left_join(meta_combined %>% dplyr::select(formula, HC, OC)) %>% 
  filter(! horizon %in% "B") %>% 
  filter(n == 1) %>% 
  gg_vankrev(aes(x = OC, y = HC, color = transect))+
  facet_wrap(~region + transect)


###############
###############

relabund_cores_neg = compute_relabund_cores(data_longform = data_negative_longform, meta_combined)
relabund_cores_pos = compute_relabund_cores(data_longform = data_positive_longform, meta_combined)

plot_relabund(relabund_cores_neg, TREATMENTS = quos(sample_label, region, transect, horizon))
plot_relabund(relabund_cores_pos, TREATMENTS = quos(sample_label, region, transect, horizon))

