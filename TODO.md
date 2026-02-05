# Overall To-do list

- Standard derivatives structures of Options as linear combination of European call puts
- Convertible pricing using more than BS
- American Call/puts
- Asian Call/puts
- All around multi-threaded Monte Carlo -> adapt to American and Asian options
- Auto-callables pricing
- Greeks calculation for exotics based on Monte-carlo
- API for live info of option prices and testing compared to market data on certain SPY options
- Web scrapping/api access for: Automatic Stock data info,Vol calculation, and Risk-free rate pull


## Front-end to do

* Ticker research, before choosing between single option characteristics input, or a way to change to exotics or option spreads and input the things needed basd on the structure with automatic pull of info from API, except date and strike. (have a toggle forallowing manual input), And after that a long or short toggle, default being long (top part of the website)
* Comparison of the pricing methods outputs in a clean table with standard divations where possible , with either my volatility calculations or the implied volatility (you can choose which one default IV) (when scrolling down)
* Vizualisation of the payout by stock price, considering the cost of the options that just was calculated (continuing to scroll down)
* Value of the Greeks for that option/position (side by side with the visualization)
* Visualization of a 1% sample of the monte-carlo paths (when continuing to scroll down)
* Market statistics section and visualization with market statistics (price, vol, price history), IV smile, and selection of option prices available through the alpha vantage API and Finnhub (at the bottom)
