

dfall = read.csv("logBFsandKSvaluesforshortvshort.csv")

plot1 = ggplot(data = dfall, mapping = aes(p, CVBFAvg)) +
  #  geom_point() +
  # geom_point(aes(p, PTC), color = "purple") +
  # geom_point(aes(p, PTN), color = "green") +
  # geom_point(aes(p, logKSb), color = "red") +
  geom_smooth(aes(p, PTC), color = "pink", se = FALSE) +
  geom_smooth(aes(p, PTN), color = "yellow", se = FALSE) +
  geom_smooth(aes(p, CVBFAvg), color = "blue", se = FALSE) +
  geom_smooth(aes(p, logKSb), color = "orange", se = FALSE) +
  labs(x = "p", y = "log(BF)") + ylim(c(-25, 45))

plot1
ggsave("BayesSimPlots/jointlogBFShortvShort.pdf",plot = plot1, device = "pdf")

plot2 = ggplot(data = dfall, mapping = aes(p, CVBFAvg)) +
  geom_point() +
  geom_point(aes(p, PTC), color = "purple") +
  geom_point(aes(p, PTN), color = "green") +
  geom_point(aes(p, logKSb), color = "red") +
  geom_smooth(aes(p, PTC), color = "pink", se = FALSE) +
  geom_smooth(aes(p, PTN), color = "yellow", se = FALSE) +
  geom_smooth(aes(p, CVBFAvg), color = "blue", se = FALSE) +
  geom_smooth(aes(p, logKSb), color = "orange", se = FALSE) +
  labs(x = "p", y = "log(BF)") + ylim(c(-70, 50))

plot2
ggsave("BayesSimPlots/jointlogBFLongvLongfull.pdf",plot = plot2, device = "pdf")