pacman::p_load(DataExplorer, gam, lubridate, GGally, MASS, tidyverse, patchwork)
set.seed(42)


tema.graficos <- theme_light() + 
  theme(plot.margin=unit(rep(.5, 4),"cm"),
        axis.title.y=element_text(margin=unit(rep(.5, 4), "cm"), size=20),
        axis.title.x=element_text(margin=unit(rep(.5, 4), "cm"), size=20),
        axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        panel.border=element_rect(colour=gray(.5)),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        legend.margin=margin(l=20, unit='pt'),
        strip.text.x=element_text(size=18),
        strip.text.y=element_text(size=18)
  )
theme_set(tema.graficos)


dados_original <- read.csv("~/data/brtvrad.csv", 
                    encoding="UTF-8",
                  stringsAsFactors=T) %>% 
  as_tibble()

dados <- dados_original %>% 
  filter(metodo != "" & 
           qt_segundo_midia > 0 &
           !is.na(transcricao_referencia) &
           transcricao_referencia != "",
           !is.na(transcricao_metodo)) %>% 
  separate(metodo, 
           into=c("metodo", "filtro"),
           sep="_") %>% 
  replace_na(list(filtro = "sem")) %>% 
  filter((metodo %in% c("audimus", "w2v") & filtro=="sem") |  
           (metodo %in% c("kaldi", "google", "azure") & filtro == "nr")) %>% 
  dplyr::select(!c(tp_gravador, cd_formato, mean_volume_nr:snr_rnnoise)) %>% 
  mutate(across(where(is.character), as.factor),
         dta_not = as_date(dta_not),
         words = n_palavras_referencia,
         # erros = pmin(n_erros_metodo, words),
         erros = pmin(n_erros_metodo, words),
         acertos = words - erros,
         wer = erros/words,
         sg_uf = fct_lump_n(sg_uf, 2, other_level="Outros"),
         filtro = fct_relevel(filtro, "sem"),
         metodo = fct_relevel(metodo, "azure")) %>%  
  rowwise() %>% 
  mutate(bitrate_kbs_factor = ifelse(bitrate_kbs > 100, "High", "Low"),
         snr_trimmed = ifelse(snr > .05, .05, snr), 
         snr_trimmed = ifelse(snr_trimmed < -.05, -.05, snr_trimmed),
         max_volume_factor = ifelse(max_volume > -2, "High", "Low")) %>% 
  ungroup() %>% 
  dplyr::select(!c(n_palavras_referencia, n_erros_metodo, hash))


mod2 <- glm(cbind(acertos, erros) ~
              (metodo +
                 sg_uf +
                 tp_veiculo +
                 bitrate_kbs_factor +
                 qt_segundo_midia +
                 mean_volume +
                 max_volume +
                 snr)^2 +
              I(qt_segundo_midia^2) +
              I(mean_volume^2) +
              I(max_volume^2) +
              I(snr^2),
      data = dados,
      family = binomial(link = "logit"),
      weights = words)

mod <- stepAIC(mod2, 
               scope = list(upper = ~
                              (metodo +
                                 sg_uf +
                                 tp_veiculo +
                                 bitrate_kbs_factor +
                                 qt_segundo_midia +
                                 mean_volume +
                                 max_volume +
                                 snr)^2 +
                              I(qt_segundo_midia^2) +
                              I(mean_volume^2) +
                              I(max_volume^2) +
                              I(snr^2),
                            lower = ~ 1),
               k=log(nrow(dados)),
               direction="backward",
               trace=T)


# dados_aux <- dados %>%
#   mutate(previsto = 1 - predict(mod, ., type="response"),
#          metodo=fct_relevel(metodo, "google", "kaldi", "audimus", "w2v", "azure"))
# levels(dados_aux$metodo) <- str_to_title(levels(dados_aux$metodo))
# 

dados_aux <- dados %>% 
  mutate(previsto = 1 - predict(mod, ., type = c("response")),
         metodo = fct_relevel(metodo, "kaldi", "google", "audimus", "w2v", "azure")) 
levels(dados_aux$metodo) <- str_to_title(levels(dados_aux$metodo))

figuras <- list(Performance=NA, Comparison=NA)

figuras[[1]] <- dados_aux %>% 
  arrange(words) %>% 
  ggplot(aes(x=previsto, 
             y=wer, 
             size=words, 
             alpha=words)) +
  geom_point(aes(color=metodo)) +
  ylab("Observed WER") + 
  xlab("Expected WER") +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by=.2), 
                     expand=expansion(mult = .02)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by=.2), 
                     expand=expansion(mult = .02)) +
  scale_color_manual(values = c("#cc3300", "#ff9966", "#ffcc00", "#99cc33", "#339900"),
                     name = "Method") +
  scale_alpha_continuous(name = "Words") +
  scale_size_continuous(name = "Words") +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  coord_fixed() + 
  geom_smooth(color = "black", show.legend=F)

# ggsave("Model_performance.pdf", width = 8, height = 8)


dados.grafico <- dados_aux %>%
  dplyr::select(cd_noticia, previsto, metodo) %>%
  pivot_wider(names_from = metodo,
              values_from = previsto) %>%
  relocate(cd_noticia, Kaldi, Google, Audimus, W2v, Azure)

medias <- dados.grafico %>% 
  colMeans(., na.rm = T)

figuras[[2]] <- dados.grafico %>%
  mutate(cd_noticia = factor(cd_noticia)) %>% 
  ggparcoord(columns = 2:6,
             groupColumn = 1,
             scale = "globalminmax",
             alphaLines = .1) +
  ylim(c(0, 1)) +
  xlab("") +
  ylab("Expected WER") + 
  scale_x_discrete(expand=expansion(mult = .05)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by=.2), 
                     expand=expansion(mult = .0)) +
  scale_color_manual(values=rep(1, length(unique(dados_aux$cd_noticia)))) +
  theme(legend.position="none") +
  geom_segment(x=1, y=medias[2], xend=2, yend=medias[3], color="#019393", size=1.3, lineend="round") +
  geom_segment(x=2, y=medias[3], xend=3, yend=medias[4], color="#019393", size=1.3, lineend="round") +
  geom_segment(x=3, y=medias[4], xend=4, yend=medias[5], color="#019393", size=1.3, lineend="round") +
  geom_segment(x=4, y=medias[5], xend=5, yend=medias[6], color="#019393", size=1.3, lineend="round") 

# ggsave("Methods_comparison.pdf", width = 8, height = 8)


figuras[[1]] + figuras[[2]]
ggsave("Figures.pdf", width = 16, height = 8)

medias_tibble  <- dados_aux %>% 
  group_by(metodo) %>% 
  summarize(wer_medio = mean(previsto)) %>% 
  arrange(wer_medio)


dados_aux %>% 
  filter(metodo %in% c("Audimus", "W2v")) %>% 
  select(cd_noticia, metodo, previsto) %>% 
  pivot_wider(names_from = metodo,
              values_from = previsto) %>% 
  rowwise() %>% 
  mutate(Audimus_best = Audimus < W2v) %>% 
  ungroup() %>% 
  pull(Audimus_best) %>% 
  mean(., na.rm=T)


dados_aux %>% 
  filter(metodo %in% c("Kaldi", "Google")) %>% 
  select(cd_noticia, metodo, previsto) %>% 
  pivot_wider(names_from = metodo,
              values_from = previsto) %>% 
  rowwise() %>% 
  mutate(Google_best = Google < Kaldi) %>% 
  ungroup() %>% 
  pull(Google_best) %>% 
  mean(., na.rm=T)
  
