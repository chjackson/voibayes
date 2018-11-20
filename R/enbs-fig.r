## Plot to illustrate Expected Net Benefit of Sampling, in manuscript.

source("enbs.r")

library(gridExtra)

p1 <- p +
  geom_line(aes(y=ec), col="red") +
  geom_line(aes(y=eb), col="blue") +
  geom_line(aes(y=eb*2), col=col2) +
  geom_line(aes(y=eb/2), col=col2) +
  geom_text(aes(x=200, y=8500), label="Benefit", col="blue", size=5) + 
  geom_text(aes(x=300, y=15000), label="Benefit (willing to pay twice)", col=col2, size=5) + 
  geom_text(aes(x=200, y=2000), label="Benefit (willing to pay half)", col=col2, size=5) + 
  geom_text(aes(x=330, y=6500), label="Cost", col="red", size=5) + 
  xlab("") +
  ylab("Cost or benefit of sampling (£)")

p2 <- p +
  geom_line(aes(y=enbs), col="black") +
  geom_line(aes(y=eb*2 - ec), col="gray") +
  geom_line(aes(y=eb/2 - ec), col="gray") +
  xlab("Additional sample size") +
  ylab("Expected net benefit of sampling (£)") +
  geom_text(aes(x=300, y=-2000), label="Willing to pay half", col="gray", size=5) + 
  geom_text(aes(x=300, y=13000), label="Willing to pay twice", col="gray", size=5) + 
  geom_vline(aes(xintercept=maxbase), col="black", linetype=2) + 
  geom_vline(aes(xintercept=maxtwice), col="gray", linetype=2) +
  geom_vline(aes(xintercept=maxhalf), col="gray", linetype=2)


# pdf("enbs.pdf")

grid.arrange(p1, p2, ncol=1)
